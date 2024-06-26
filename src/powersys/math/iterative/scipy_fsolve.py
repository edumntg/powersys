from .iterative import *
from ...models import PowerSystem
import numpy as np
from scipy.optimize import fsolve

class ScipyFsolve(Iterative):
    def __init__(self, model: PowerSystem, args: IterativeArgs):
        super().__init__(model, args)

    def solve(self, state_dict, disp = False):
        # Create initial guess
        nvars = 2*self.model.n_buses

        # Create vector of initial values
        x0 = np.zeros(nvars)
        v = 0
        for bus in self.model.buses:
            if bus.type == PowerSystem.SLACK: # variables are P and Q
                x0[v] = bus.Pgen
                x0[v+1] = bus.Qgen
            elif bus.type == PowerSystem.PV: # variables are angle and Qgen
                x0[v] = bus.angle
                x0[v+1] = bus.Qgen
            else: # bus is PQ
                x0[v] = bus.V
                x0[v+1] = bus.angle

            v += 2
        
        # Solve
        x = fsolve(self.__scipy_func, x0, xtol = self.tol)
        
        # Re-arrange values into vectors
        V = state_dict['variables']['V']
        theta =state_dict['variables']['theta']
        Pgen = state_dict['variables']['Pgen']
        Qgen = state_dict['variables']['Qgen']

        v = 0
        for bus in self.model.buses:
            if bus.type == PowerSystem.SLACK: # slack: vars are P and Q
                Pgen[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            elif bus.type == PowerSystem.PV: # pv: vars are theta and Q
                theta[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            else: # pq: vars are V and theta
                V[bus.id] = x[v]
                theta[bus.id] = x[v+1]
        
            v += 2

        return V, theta, Pgen, Qgen

    def __scipy_func(self, x):
        V = np.array([bus.V for bus in self.model.buses])
        theta = np.array([bus.angle for bus in self.model.buses])
        Pgen = np.array([bus.Pgen for bus in self.model.buses])
        Qgen = np.array([bus.Qgen for bus in self.model.buses])

        v = 0
        for bus in self.model.buses:
            if bus.type == PowerSystem.SLACK: # slack: vars are P and Q
                Pgen[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            elif bus.type == PowerSystem.PV: # pv: vars are theta and Q
                theta[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            else: # pq: vars are V and theta
                V[bus.id] = x[v]
                theta[bus.id] = x[v+1]
        
            v += 2
        
        Pout = np.zeros(self.model.n_buses)
        Qout = np.zeros(self.model.n_buses)

        for line in self.model.lines:
            i = line.from_bus
            k = line.to_bus
            if i == k:
                continue

            Pout[i] += ((-self.model.G[i,k] + self.model.g[i,k])*V[i]**2 + V[i]*V[k]*(self.model.G[i,k]*np.cos(theta[i] - theta[k]) + self.model.B[i,k]*np.sin(theta[i] - theta[k])))
            Pout[k] += ((-self.model.G[k,i] + self.model.g[k,i])*V[k]**2 + V[k]*V[i]*(self.model.G[k,i]*np.cos(theta[k] - theta[i]) + self.model.B[k,i]*np.sin(theta[k] - theta[i])))

            Qout[i] += ((self.model.B[i,k] - self.model.b[i,k])*V[i]**2 + V[i]*V[k]*(-self.model.B[i,k]*np.cos(theta[i] - theta[k]) + self.model.G[i,k]*np.sin(theta[i] - theta[k])))
            Qout[k] += ((self.model.B[k,i] - self.model.b[k,i])*V[k]**2 + V[k]*V[i]*(-self.model.B[k,i]*np.cos(theta[k] - theta[i]) + self.model.G[k,i]*np.sin(theta[k] - theta[i])))

        # Finally, set equations
        equations = []
        for bus in self.model.buses:

            equations.append(
                Pgen[bus.id] + bus.Pgen_fixed - bus.Pload - Pout[bus.id]
            )
            equations.append(
                Qgen[bus.id] + bus.Qgen_fixed - bus.Qload - Qout[bus.id]
            )

        return equations