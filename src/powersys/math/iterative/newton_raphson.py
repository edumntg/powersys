from .iterative import *
from ...models import PowerSystem
import numpy as np
from tabulate import tabulate

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

        results_table = []
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
            m = 0
            for i, bus_i in enumerate(self.model.buses): # This loop will contain all theta variables with respect with are calculating derivatives
                if bus_i.type == PowerSystem.SLACK: # ignore slack bus because theta changes only for PV and PQ buses
                    continue
                
                n = 0
                for k, bus_k in enumerate(self.model.buses): # This loop controls the variable with respect to which we are calculating the derivative (theta2, theta3...)
                    # We are calculating the derivative with respect to buses that have variable angles. That is, we ignore slack buses
                    if bus_k.type == PowerSystem.SLACK:
                        continue
                        
                    # Now, we use an additional loop to calculate all flows from i to buses 1..n
                    if i == k: # All terms of Pi are dependent of theta_k, so all derivatives are ~= 0
                        for j, bus_j in enumerate(self.model.buses):
                            J1[m,n] += V[i]*V[j]*(-G[i,j]*np.sin(theta[i] - theta[j]) + B[i,j]*np.cos(theta[i] - theta[j]))

                        J1[m,n] += -V[i]**2 * B[i,i]

                    else: # Only term for j == k is different from zero
                        J1[m,n] = V[i]*V[j]*(G[i,j]*np.sin(theta[i] - theta[j]) - B[i,j]*np.cos(theta[i] - theta[j]))
                    
                    n += 1
                
                m += 1

            # J2 will contain all derivatives of P with respect to voltages
            m = 0
            for i, bus_i in enumerate(self.model.buses): # This loop will contain all voltage variables with respect with are calculating derivatives
                if bus_i.type == PowerSystem.SLACK:
                    continue
                
                n = 0
                for k, bus_k in enumerate(self.model.buses): # This loop controls the variable with respect to which we are calculating the derivative (V2, V3...)
                    # We are calculating the derivative with respect to buses that have variable voltages. That is, we ignore slack/pv buses
                    if bus_k.type != PowerSystem.PQ:
                        continue
                        
                    # Now, we use an additional loop to calculate all flows from i to buses 1..n
                    if i == k: # All terms of Pi are dependent of theta_k, so all derivatives are ~= 0
                        for j, bus_j in enumerate(self.model.buses):
                            J2[m,n] += 2*V[i]*(-G[i,j] + g[i,j]) + V[j]*(G[i,j]*np.cos(theta[i] - theta[j]) + B[i,j]*np.sin(theta[i] - theta[j]))
                        
                        J2[m,n] += V[i]*G[i,i]

                    else: # Only term for j == k is different from zero
                        J2[m,n] = V[i]*(G[i,j]*np.cos(theta[i] - theta[j]) + B[i,j]*np.sin(theta[i] - theta[j]))
                    
                    n += 1
                
                m += 1

            # J3 will contain all derivatives of Q with respect to angles
            m = 0
            for i, bus_i in enumerate(self.model.buses): # This loop will contain all theta variables with respect with are calculating derivatives
                if bus_i.type != PowerSystem.PQ: # ignore slack/pv buses
                    continue
                
                n = 0
                for k, bus_k in enumerate(self.model.buses): # This loop controls the variable with respect to which we are calculating the derivative (theta2, theta3...)
                    # We are calculating the derivative with respect to buses that have variable angles. That is, we ignore slack buses
                    if bus_k.type == PowerSystem.SLACK:
                        continue
                        
                    # Now, we use an additional loop to calculate all flows from i to buses 1..n
                    if i == k: # All terms of Pi are dependent of theta_k, so all derivatives are ~= 0
                        for j, bus_j in enumerate(self.model.buses):
                            J3[m,n] += V[i]*V[j]*(B[i,j]*np.sin(theta[i] - theta[j]) + G[i,j]*np.cos(theta[i] - theta[j]))

                        J3[m,n] += -V[i]**2 *G[i,i]

                    else: # Only term for j == k is different from zero
                        J3[m,n] = V[i]*V[j]*(-B[i,j]*np.sin(theta[i] - theta[j]) - G[i,j]*np.cos(theta[i] - theta[j]))
                    
                    n += 1

                m += 1

            # J4 will contain all derivatives of Q with respect to voltages
            m = 0
            for i, bus_i in enumerate(self.model.buses): # This loop will contain all theta variables with respect with are calculating derivatives
                if bus_i.type != PowerSystem.PQ: # ignore slack/pv buses
                    continue
                
                n = 0
                for k, bus_k in enumerate(self.model.buses): # This loop controls the variable with respect to which we are calculating the derivative (V2, V3...)
                    # We are calculating the derivative with respect to buses that have variable voltages. That is, we ignore slack/pv buses
                    if bus_k.type != PowerSystem.PQ:
                        continue
                        
                    # Now, we use an additional loop to calculate all flows from i to buses 1..n
                    if i == k: # All terms of Pi are dependent of theta_k, so all derivatives are ~= 0
                        for j, bus_j in enumerate(self.model.buses):
                            J4[m,n] += 2*V[i]*(B[i,j] - b[i,j]) + V[j]*(-B[i,j]*np.cos(theta[i] - theta[j]) + G[i,j]*np.sin(theta[i] - theta[j]))

                        J4[m,n] += -V[i]*B[i,i]

                    else: # Only term for j == k is different from zero
                        J4[m,n] = V[i]*(-B[i,j]*np.cos(theta[i] - theta[j]) + G[i,j]*np.sin(theta[i] - theta[j]))
                    
                    n += 1
                
                m += 1

            # Construct Jacobian
            J = np.vstack((
                np.hstack((J1, J2)),
                np.hstack((J3, J4))
            ))

            # Compute mismatches
            deltaP = Psp - P
            deltaQ = Qsp - Q

            # Vector of mismatches contains only mismatches for PV and PQ variables
            dPQ = np.zeros((self.model.n_pv + 2*self.model.n_pq, 1))
            m = 0
            for i, bus in enumerate(self.model.buses):
                if bus.type == PowerSystem.SLACK:
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

            dtheta = dx[:self.model.N-1]
            dV = dx[self.model.N-1:]

            
            # Update angle for all non-slack buses
            x[1:self.model.N] += dtheta

            # Update voltage for all PQ buses
            x[-self.model.n_pq:] += dV

            # Compute errors
            err = np.max(np.abs(dx))

            results_table.append([iters, err])

        if self.args.verbose:
            print(tabulate(results_table, headers = ["Iter No.", "Error"]))

        if iters < self.max_iters:
            print(f"Load flow solved in {iters} iterations!")
        else:
            print(f"Max. number of iterations reached. Load flow did not converge!")
            raise Exception("Load Flow did not converge!")

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

        return np.abs(V)[:,0], np.angle(V)[:,0], P[:,0], Q[:,0]