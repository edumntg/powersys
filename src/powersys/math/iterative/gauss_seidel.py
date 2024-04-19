from .iterative import *
from ...models.powersystem import *
import numpy as np

class GaussSeidel(Iterative):
    def __init__(self, args: IterativeArgs):
        super().__init__(args)

    def solve(self, model: PowerSystem):
        V = np.array([
            np.abs(bus.V)*(np.cos(bus.theta) + 1j*np.sin(bus.theta))
            for bus in model.buses
        ], dtype="complex")

        P = np.array([bus.Pgen - bus.Pload for bus in model.buses])
        Q = np.array([bus.Qgen - bus.Qload for bus in model.buses])

        iters = 1
        err = 1E9
        Vprev = V.copy()
        while err > self.tol and iters < self.max_iters:
            for bus in model.buses:
                if bus.reference: # slack
                    continue

                PYV = 0
                for bus2 in model.buses:
                    if bus.id != bus2.id:
                        PYV += model.Ybus[bus.id, bus2.id]*V[bus.id]
                
                is_pv = model.generators.some(lambda x: x.bus == bus.id)
                if is_pv:
                    Q[bus.id] = -np.imag(np.conj(V[bus.id])*(PYV + model.Ybus[bus.id, bus.id]*V[bus.id]))

                # Compute new voltages
                V[bus.id] = (1/model.Ybus[bus.id, bus.id])*((P[bus.id] - 1j*Q[bus.id])/np.conj(V[bus.id]) - PYV)

                if is_pv:
                    # Maintain voltage amplite
                    V[bus.id] = np.abs(Vprev[bus.id])*(np.cos(np.angle(V[bus.id])) + 1j*np.sin(np.angle(V[bus.id])))

            # Compute error
            err = np.max(np.abs(np.abs(V) - np.abs(Vprev)))
            iters += 1

            if iters >= self.max_iters:
                print(f"Max. number of iterations ({iters}) reached. Gauss-Seidel method did not converge.")
            
            Vprev = V.copy()
            