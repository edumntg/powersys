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

        # J1 will have size: (N-1, N-1)
        J1 = np.zeros((self.model.N-1, self.model.N-1))

        # J2 will have size: (N-1, n_pq)
        J2 = np.zeros((self.model.N-1, self.model.n_pq))

        # J3 will have size: (n_pq, N-1)
        J3 = np.zeros((self.model.n_pq, self.model.N-1))

        # J4 will have size: (n_pq, n_pq)
        J4 = np.zeros((self.model.n_pq, self.model.n_pq))

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

            # Reset gradients
            J1.fill(0)
            J2.fill(0)
            J3.fill(0)
            J4.fill(0)

            # J1 will contain all derivatives of P with respect to angles
            for m in range(self.model.N-1):
                n = 0
                for i, bus_i in enumerate(self.model.buses):
                    if bus_i.type == PowerSystem.SLACK: # ignore slack bus
                        continue

                    dPi_dthetak = 0
                    for k, bus_k in enumerate(self.model.buses):
                        # We compute dPdtheta with respect to angles of PV and PQ buses
                        if bus_k.type == PowerSystem.SLACK:
                            continue

                        if i == k:
                            dPi_dthetak = V[i]*V[k]*(-G[i,k]*np.sin(theta[i] - theta[k]) + B[i,k]*np.cos(theta[i] - theta[k]))
                        else:
                            dPi_dthetak = V[i]*V[k]*(G[i,k]*np.sin(theta[i] - theta[k]) - B[i,k]*np.cos(theta[i] - theta[k]))

                    J1[m,n] = dPi_dthetak
                    n += 1

            # J2 will contain all derivatives of P with respect to voltages
            for m in range(self.model.N-1):
                n = 0
                for i, bus_i in enumerate(self.model.buses):
                    if bus_i.type != PowerSystem.PQ: # consider PQ buses only
                        continue

                    dPi_dVk = 0
                    for k, bus_k in enumerate(self.model.buses):
                        # We compute dPdV with respect to voltages of PQ buses
                        if bus_k.type != PowerSystem.PQ:
                            continue

                        if i == k:
                            dPi_dVk = 2*V[i]*(-G[i,k] + g[i,k]) + V[k]*(G[i,k]*np.cos(theta[i] - theta[k]) + B[i,k]*np.sin(theta[i] - theta[k]))
                        else:
                            dPi_dVk = V[i]*(G[i,k]*np.cos(theta[i] - theta[k]) + B[i,k]*np.sin(theta[i] - theta[k]))

                    J2[m,n] = dPi_dVk
                    n += 1

            # J3 will contain all derivatives of Q with respect to angles
            for m in range(self.model.n_pq):
                n = 0
                for i, bus_i in enumerate(self.model.buses):
                    if bus_i.type == PowerSystem.SLACK: # consider PQ buses only
                        continue

                    dQi_dthetak = 0
                    for k, bus_k in enumerate(self.model.buses):
                        # We compute dPdV with respect to voltages of PQ buses
                        if bus_k.type == PowerSystem.SLACK:
                            continue

                        if i == k:
                            dQi_dthetak = V[i]*V[k]*(B[i,k]*np.cos(theta[i] - theta[k]) + G[i,k]*np.sin(theta[i] - theta[k]))
                        else:
                            dQi_dthetak = V[i]*V[k]*(-B[i,k]*np.sin(theta[i] - theta[k]) - G[i,k]*np.cos(theta[i] - theta[k]))

                    J3[m,n] = dQi_dthetak
                    n += 1

            # J4 will contain all derivatives of Q with respect to voltages
            for m in range(self.model.n_pq):
                n = 0
                for i, bus_i in enumerate(self.model.buses):
                    if bus_i.type != PowerSystem.PQ: # consider PQ buses only
                        continue

                    dQi_dVk = 0
                    for k, bus_k in enumerate(self.model.buses):
                        # We compute dQdV with respect to voltages of PQ buses
                        if bus_k.type != PowerSystem.PQ:
                            continue

                        if i == k:
                            dQi_dVk = 2*V[i]*(B[i,k] - b[i,k]) + V[k]*(-B[i,k]*np.cos(theta[i] - theta[k]) + G[i,k]*np.sin(theta[i] - theta[k]))
                        else:
                            dQi_dVk = V[i]*(-B[i,k]*np.cos(theta[i] - theta[k]) + G[i,k]*np.sin(theta[i] - theta[k]))

                    J4[m,n] = dQi_dVk
                    n += 1

            # Construct Jacobian
            J = np.vstack((
                np.hstack((J1, J2)),
                np.hstack((J3, J4))
            ))
            print(J)

            # Compute mismatches
            deltaP = Psp - P
            deltaQ = Qsp - Q

            # Vector of mismatches contains only mismatches for PV and PQ variables
            dPQ = np.zeros((self.model.n_pv + 2*self.model.n_pq, 1))
            m = 0
            for i, bus in enumerate(self.model.buses):
                if bus.type != PowerSystem.PV:
                    continue

                dPQ[m] = deltaP[i]
                m += 1

            for i, bus in enumerate(self.model.buses):
                if bus.type != PowerSystem.PQ:
                    continue

                dPQ[m] = deltaQ[i]
                m += 1
        
            # Compute corrections in variables
            dx = np.linalg.solve(J, dPQ) # [dtheta, dV]

            # Compute new variables
            m = 0
            # Update angles
            for i, bus in enumerate(self.model.buses):
                if bus.type != PowerSystem.SLACK:
                    x[i] += dx[m]
                    m += 1

            # Update voltages
            for i, bus in enumerate(self.model.buses):
                if bus.type == PowerSystem.PQ:
                    x[i] += dx[m]
                    m += 1 

            # Compute errors
            err = np.max(np.abs(dx))

            # For slack and PV buses, fix voltages
            # for i, bus in enumerate(self.model.buses):
            #     if bus.type == PowerSystem.SLACK:
            #         x[self.model.N+i] = bus.V
            #         x[i] = bus.angle
            #     elif bus.type == PowerSystem.PV:
            #         x[self.model.N+i] = bus.V
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