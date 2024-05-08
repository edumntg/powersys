from .iterative import *
from ...models import PowerSystem

import numpy as np

class NewtonRaphson(Iterative):
    def __init__(self, model: PowerSystem, args: IterativeArgs):
        super().__init__(model, args)

    def solve(self, state_dict, disp = False):
        # For Newton-Raphson, the variables are Voltages and angles
        variables = state_dict['variables']
        Vmag = variables['V']
        theta = variables['theta']

        # We create a vector of size 2*n containing the initial values of the variables ( theta1, theta2, ... thetan, v1, v2,.. vn)
        x = np.concatenate((np.array([bus.theta for bus in self.model.buses]), np.array([bus.V for bus in self.model.buses]))).reshape((2*self.model.N, 1))

        iters = 0
        err = 1E9

        G = self.model.G
        B = self.model.B
        g = self.model.g
        b = self.model.b


        while err > self.tol and iters < self.max_iters:
            iters += 1

            # Get current voltages in complex form
            theta = x[:self.model.N]
            V = x[self.model.N:]

            # Compute powers
            P = np.zeros((self.model.N, 1))
            Q = np.zeros((self.model.N, 1))
            for i, bus_i in enumerate(self.model.buses):
                for k, bus_k in enumerate(self.model.buses):
                    P[i] += (-G[i,k] + g[i,k])*V[i]**2 + V[i]*V[k]*(G[i,k]*np.cos(theta[i] - theta[k]) + B[i,k]*np.sin(theta[i] - theta[k]))
                    Q[i] += (B[i,k] - b[i,k])*V[i]**2 + V[i]*V[k]*(-B[i,k]*np.cos(theta[i] - theta[k]) + G[i,k]*np.sin(theta[i] - theta[k]))

            # Compute setpoints
            Psp = np.array([bus.Pgen - bus.Pload for bus in self.model.buses]).reshape((self.model.N, 1))
            Qsp = np.array([bus.Qgen - bus.Qload for bus in self.model.buses]).reshape((self.model.N, 1))

            # Create terms of the Jacobian

            # J1 contains derivatives of P with respect to angles
            J1 = np.zeros((self.model.N, self.model.N))
            # J2 contains derivatives of P with respect to voltages
            J2 = np.zeros((self.model.N, self.model.N))
            # J3 contains derivatives of Q with respect to angles
            J3 = np.zeros((self.model.N, self.model.N))
            # J4 contains derivatives of Q with respect to voltages 
            J4 = np.zeros((self.model.N, self.model.N))

            # J1 contains the derivatives of P with respect to angles
            # J2 contains the derivatives of P with respect to voltages
            # J3 contains the derivatives of Q with respect to angles
            # J4 contains the derivatives of Q with respect to voltages
            for i, bus_i in enumerate(self.model.buses):
                for k, bus_k in enumerate(self.model.buses):
                    if i == k:
                        #dP_dthetai
                        J1[i,k] += V[i]*V[k]*(-G[i,k]*np.sin(theta[i]-theta[k]) + B[i,k]*np.cos(theta[i]-theta[k]))
                        #dP_dVi
                        J2[i,k] += 2*V[i]*(-G[i,k] + g[i,k]) + V[k]*(G[i,k]*np.cos(theta[i]-theta[k]) + B[i,k]*np.sin(theta[i]-theta[k]))
                        #dQ_dthetai
                        J3[i,k] += V[i]*V[k]*(B[i,k]*np.sin(theta[i]-theta[k]) + G[i,k]*np.cos(theta[i]-theta[k]))
                        #dQ_dVi
                        J4[i,k] += 2*V[i]*(B[i,k]-b[i,k]) + V[k]*(-B[i,k]*np.cos(theta[i]-theta[k]) + G[i,k]*np.sin(theta[i]-theta[k]))
                    else:
                        #dP_dthetak
                        J1[i,k] += V[i]*V[k]*(G[i,k]*np.sin(theta[i]-theta[k]) - B[i,k]*np.cos(theta[i]-theta[k]))
                        #dP_dVk
                        J2[i,k] += V[i]*(G[i,k]*np.cos(theta[i]-theta[k]) + B[i,k]*np.sin(theta[i]-theta[k]))
                        #dQ_dthetak
                        J3[i,k] += V[i]*V[k]*(-B[i,k]*np.sin(theta[i]-theta[k])-G[i,k]*np.cos(theta[i]-theta[k]))
                        #dQ_dVk
                        J4[i,k] += V[i]*(-B[i,k]*np.cos(theta[i]-theta[k]) + G[i,k]*np.sin(theta[i]-theta[k]))

            # Form J
            J = np.vstack((np.hstack((J1, J2)), np.hstack((J3, J4))))

            # Compute mismatches
            deltaP = Psp - P
            deltaQ = Qsp - Q
            
            deltaf = np.vstack((deltaP, deltaQ))

            # Compute corrections in variables
            dx = np.linalg.solve(J, deltaf) # [dtheta, dV]

            # Compute new variables
            x_new = x + dx

            # Compute errors
            err = np.max(np.abs(dx))

            # Update current values
            x = x_new.copy()

            # For slack and PV buses, fix voltages
            for i, bus in enumerate(self.model.buses):
                if bus.type == PowerSystem.SLACK:
                    x[self.model.N+i] = bus.V
                    x[i] = bus.angle
                elif bus.type == PowerSystem.PV:
                    x[self.model.N+i] = bus.V
            # Print error
            print(f"Iteration {iters}, error {err}")

        theta = x[:self.model.N]
        V = x[self.model.N:]

        V = V*(np.cos(theta) + 1j*np.sin(theta))
        I = self.model.Ybus@V
        S = V*np.conj(I)
        P = np.real(S)
        Q = np.imag(S)

        for bus in self.model.buses:
            if bus.type == PowerSystem.SLACK: # slack
                P[bus.id] = np.real(S[bus.id])
                Q[bus.id] = np.imag(S[bus.id])
            elif bus.type == PowerSystem.PV: #pv
                P[bus.id] = np.real(S[bus.id]) + bus.Pload
                Q[bus.id] = np.imag(S[bus.id]) + bus.Qload

        return np.abs(V), np.angle(V), P, Q