from .iterative import *
from ...models import PowerSystem
import numpy as np
from scipy.optimize import fsolve

class ScipyFsolve(Iterative):
    def __init__(self, model: PowerSystem, args: IterativeArgs):
        super().__init__(model, args)

    def solve(self, disp = False):
        # Create initial guess
        nvars = 2*self.model.n_buses
        x0 = np.random.randn(nvars)
        
        # Solve
        x = fsolve(self.__scipy_func, x0)
        
        # Re-arrange values into vectors
        V = np.array([bus.V for bus in self.model.buses])
        theta = np.array([bus.angle for bus in self.model.buses])
        Pgen = np.array([bus.Pgen for bus in self.model.buses])
        Qgen = np.array([bus.Qgen for bus in self.model.buses])

        v = 0
        for bus in self.model.buses:
            if bus.type == 1: # slack: vars are P and Q
                Pgen[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            elif bus.type == 2: # pv: vars are theta and Q
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
            if bus.type == 1: # slack: vars are P and Q
                Pgen[bus.id] = x[v]
                Qgen[bus.id] = x[v+1]
            elif bus.type == 2: # pv: vars are theta and Q
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
            Pout[i] += ((-self.model.G[i,k] + self.model.g[i,k])*V[i]**2 + V[i]*V[k]*(self.model.G[i,k]*np.cos(theta[i] - theta[k]) + self.model.B[i,k]*np.sin(theta[i] - theta[k])))
            Pout[k] += ((-self.model.G[k,i] + self.model.g[k,i])*V[k]**2 + V[k]*V[i]*(self.model.G[k,i]*np.cos(theta[k] - theta[i]) + self.model.B[k,i]*np.sin(theta[k] - theta[i])))

            Qout[i] += ((self.model.B[i,k] - self.model.b[i,k])*V[i]**2 + V[i]*V[k]*(-self.model.B[i,k]*np.cos(theta[i] - theta[k]) + self.model.G[i,k]*np.sin(theta[i] - theta[k])))
            Qout[k] += ((self.model.B[k,i] - self.model.b[k,i])*V[k]**2 + V[k]*V[i]*(-self.model.B[k,i]*np.cos(theta[k] - theta[i]) + self.model.G[k,i]*np.sin(theta[k] - theta[i])))

        # Finally, set equations
        equations = []
        for bus in self.model.buses:

            equations.append(
                Pgen[bus.id] - bus.Pload - Pout[bus.id]
            )
            equations.append(
                Qgen[bus.id] - bus.Qload - Qout[bus.id]
            )

        return equations