from .iterative import *
from ...models import PowerSystem
import numpy as np
from tabulate import tabulate

class GaussSeidel(Iterative):
    def __init__(self, model: PowerSystem, args: IterativeArgs):
        super().__init__(model, args)

    def solve(self, state_dict, disp = False):
        variables = state_dict['variables']
        Vmag = variables['V']
        theta = variables['theta']
        Pgen = variables['Pgen']
        Qgen = variables['Qgen']
        Pgen_fixed = variables['Pgen_fixed']
        Qgen_fixed = variables['Qgen_fixed']
        V = np.array([
            np.abs(Vmag[bus.id])*(np.cos(theta[bus.id]) + 1j*np.sin(theta[bus.id]))
            for bus in self.model.buses
        ], dtype="complex")

        P = np.array([bus.Pgen - bus.Pload for bus in self.model.buses])
        Q = np.array([bus.Qgen - bus.Qload for bus in self.model.buses])

        iters = 1
        err = 1E9
        Vprev = V.copy()

        results_table = []
        while err > self.tol and iters < self.max_iters:
            for bus in self.model.buses:
                if bus.type == PowerSystem.SLACK: # slack
                    continue

                PYV = 0
                for bus2 in self.model.buses:
                    if bus.id != bus2.id:
                        PYV += self.model.Ybus[bus.id, bus2.id]*V[bus2.id]
                
                is_pv = self.model.generators.some(lambda x: x.bus == bus.id)
                if is_pv:
                    Q[bus.id] = -np.imag(np.conj(V[bus.id])*(PYV + self.model.Ybus[bus.id, bus.id]*V[bus.id]))

                # Compute new voltages
                V[bus.id] = (1/self.model.Ybus[bus.id, bus.id])*((P[bus.id] - 1j*Q[bus.id])/np.conj(V[bus.id]) - PYV)

                if is_pv:
                    # Maintain voltage amplite
                    V[bus.id] = np.abs(Vprev[bus.id])*(np.cos(np.angle(V[bus.id])) + 1j*np.sin(np.angle(V[bus.id])))

            # Compute error
            err = np.max(np.abs(np.abs(V) - np.abs(Vprev)))
            results_table.append([iters, err])

            iters += 1

            Vprev = V.copy()

        if self.args.verbose:
            print(tabulate(results_table, headers = ["Iter No.", "Error"]))

        if iters < self.max_iters:
            print(f"Load flow solved in {iters} iterations!")
        else:
            print(f"Max. number of iterations ({iters}) reached. Gauss-Seidel method did not converge.")
            raise Exception("Load Flow did not converge!")
            
        # Compute currents
        I = self.model.Ybus@V
        S = V*np.conj(I)

        for bus in self.model.buses:
            if bus.type == PowerSystem.SLACK: # slack
                P[bus.id] = np.real(S[bus.id])
                Q[bus.id] = np.imag(S[bus.id])
            elif bus.type == PowerSystem.PV: #pv
                P[bus.id] = np.real(S[bus.id]) + bus.Pload
                Q[bus.id] = np.imag(S[bus.id]) + bus.Qload
        
        return np.abs(V), np.angle(V), P, Q
            