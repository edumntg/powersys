from ..models.PowerSystem import PowerSystem
from ..math.explicit.RK4 import RK4
import numpy as np
from gekko import GEKKO

class TransientAnalysisSolver(object):

    def __init__(self, system: PowerSystem, tspan: np.array):
        self.system = system
        self.solved = False
        self.tspan = tspan
        self.__dt = tspan[1] - tspan[0]

        self.__initialize_values() # Initialize all vectors and pre-event values

    def __initialize_values(self):
        """Initializes all vectors/arrays required to store all values during the simulation.

        Args:
            None

        Returns:
            None

        """

        if not self.system:
            raise "No PowerSystem object given"
        
        N = len(self.tspan) # Number of time points

        self.delta = np.zeros((self.system.n_gens, N))
        self.omega = np.zeros((self.system.n_gens, N))
        self.Pegap = np.zeros((self.system.n_gens, N))
        self.Eqp = np.zeros((self.system.n_gens, N))
        self.Edp = np.zeros((self.system.n_gens, N))
        self.Eqpp = np.zeros((self.system.n_gens, N))
        self.Edpp = np.zeros((self.system.n_gens, N))
        self.Vt = np.zeros((self.system.n_buses, N))
        self.theta = np.zeros((self.system.n_buses, N))
        self.Iq = np.zeros((self.system.n_gens, N))
        self.Id = np.zeros((self.system.n_gens, N))
        self.Pmgap = np.zeros((self.system.n_gens, N))
        self.Xv = np.zeros((self.system.n_gens, N))
        self.Pc = np.zeros((self.system.n_gens, N))
        self.Vmed = np.zeros((self.system.n_gens, N))
        self.Vvi = np.zeros((self.system.n_gens, N))
        self.Va = np.zeros((self.system.n_gens, N))
        self.Vexc = np.zeros((self.system.n_gens, N))
        self.Vw = np.zeros((self.system.n_gens, N))
        self.Vpc1 = np.zeros((self.system.n_gens, N))
        self.Vpc2 = np.zeros((self.system.n_gens, N))

        self.ev = np.zeros((self.system.n_gens, N))

        # Assign initial values
        self.delta[:,0] = np.array([gen.df for gen in self.system.generators])
        self.Pmgap[:,0] = np.array([gen.Pm for gen in self.system.generators])
        self.Pegap[:,0] = self.Pmgap[:,0]
        self.Iq[:,0] = np.array([gen.Iq for gen in self.system.generators])
        self.Id[:,0] = np.array([gen.Id for gen in self.system.generators])
        self.Eqp[:,0] = np.array([gen.Eqp for gen in self.system.generators])
        self.Edp[:,0] = np.array([gen.Edp for gen in self.system.generators])
        self.Eqpp[:,0] = np.array([gen.Eqpp for gen in self.system.generators])
        self.Edpp[:,0] = np.array([gen.Edpp for gen in self.system.generators])
        self.Xv[:,0] = self.Pegap[:,0]
        self.Pc[:,0] = self.Pegap[:,0]
        self.Vmed[:,0] = self.system.Vmed0
        self.Vexc[:,0] = self.system.Vexc0
        self.Vvi[:,0] = np.array([gen.Vvi for gen in self.system.generators])
        self.Va[:,0] = np.array([gen.Va for gen in self.system.generators])
        self.Vt[:,0] = np.array([bus.V for bus in self.system.buses])
        self.theta[:,0] = np.array([bus.theta for bus in self.system.generators])

        # Vectors to store diff values
        self.dw = np.zeros((self.system.n_gens, 1))
        self.dd = np.zeros((self.system.n_gens, 1))
        self.dVmed = np.zeros((self.system.n_gens, 1))
        self.ev = np.zeros((self.system.n_gens, 1))
        self.dVvi = np.zeros((self.system.n_gens, 1))
        self.dVw = np.zeros((self.system.n_gens, 1))
        self.dVpc1 = np.zeros((self.system.n_gens, 1))
        self.dVpc2 = np.zeros((self.system.n_gens, 1))
        self.dVa = np.zeros((self.system.n_gens, 1))
        self.dVexc = np.zeros((self.system.n_gens, 1))
        self.dEqp = np.zeros((self.system.n_gens, 1))
        self.dEdp = np.zeros((self.system.n_gens, 1))
        self.dEqpp = np.zeros((self.system.n_gens, 1))
        self.dEdpp = np.zeros((self.system.n_gens, 1))
        self.Xp = np.zeros((self.system.n_gens, 1))
        self.dXv = np.zeros((self.system.n_gens, 1))
        self.dPc = np.zeros((self.system.n_gens, 1))
        self.Vx = np.zeros((self.system.n_gens, 1))
        self.Vvp = np.zeros((self.system.n_gens, 1))
        self.Vvd = np.zeros((self.system.n_gens, 1))

    def __state_equations2(self, iter = 0):
        # For each generator, we will have 7 equations
        for i, gen in enumerate(self.system.generators):

            self.dw[i] = (0.5/gen.H) * (self.Pmgap[i, iter] - self.Pegap[i, iter])
            self.dd[i] = self.w[i]*self.system.we
            self.dEqp[i] = (1.0/gen.Td0p)*(-self.Eqp[i, iter] + (gen.Xd - gen.Xdp)*self.Id[i, iter] + self.Vexc[i, iter])
            self.dEdp[i] = (1.0/gen.Tq0p)*(-self.Edp[i, iter] - (gen.Xq - gen.Xqp)*self.Iq[i, iter])

            if gen.USE_CAGE:
                self.dEqpp[i] = (1.0/gen.Td0pp)*(-self.Eqpp[i, iter] + (gen.Xd - gen.Xdpp)*self.Id[i, iter] + self.Eqp[i, iter])
                self.dEdpp[i] = (1.0/gen.Tq0pp)*(-self.Edpp[i, iter] - (gen.Xq - gen.Xqpp)*self.Iq[i, iter] + self.Edp[i, iter])

            self.dVmed[i] = (1.0/gen.Tmed)*(gen.Kmed*np.abs(self.Vt[i, iter] - self.Vmed[i, iter]))
            if gen.USE_PSS:
                self.ev[i] = self.Vt[i, iter] - self.Vmed[i, iter] + self.Vpc2[i, iter]
            else:
                self.ev[i] = self.Vt[i, iter] - self.Vmed[i, iter]

            self.dVvi[i] = gen.Kv*self.ev[i]
            
            # PSS
            self.dVw[i] = gen.Kest*self.dw[i] - self.Vw[i, iter]/gen.Tw
            self.dVpc1[i] = (1/gen.T2)*(self.Vw[i, iter] + gen.T1*self.dVw[i] - self.Vpc1[i, iter])
            self.dVpc2[i] = (1/gen.T4)*(self.Vpc1[i, iter] + gen.T3*self.dVpc1[i] - self.Vpc2[i, iter])

            self.Vvp[i] = gen.Kp*self.ev[i]

            if gen.USE_PSS:
                self.Vvd[i] = gen.Kd*(-self.dVmed[i] + self.dVpc2[i])
            else:
                self.Vvd[i] = -gen.Kd*self.dVmed[i]

            self.Vx[i] = self.Vvp[i] + self.Vvd[i] + self.Vvi[i, iter]
            self.dVa[i] = (1.0/gen.Ta)*(gen.Ka*self.Vx[i] - self.Va[i, iter])
            self.dVexc[i] = (1.0/gen.Texc)*(gen.Kexc*self.Va[i, iter] - self.Vexc[i, iter])

            self.Xp[i] = -self.w[i, iter]/gen.R

            self.dXv[i] = (1.0/gen.Tv)*(gen.Kv*(self.Pc[i, iter] + gen.Xp) - self.Xv[i, iter])
            self.dPc[i] = -gen.Ki*self.w[i, iter]

    def __compute_non_dynamic_values(self, iter = 0):
        for i, gen in enumerate(self.system.generators):
            self.Eqp[i,iter] = self.Vq[i, iter] + gen.Ra*self.Iq[i, iter] - gen.Xdp*self.Id[i, iter]
            self.Edp[i,iter] = self.Vd[i, iter] + gen.Ra*self.Id[i, iter] + gen.Xqp*self.Iq[i, iter]
            self.Eqpp[i,iter] = self.Vq[i,iter] + gen.Ra*self.Iq[i, iter] - gen.Xdpp*self.Id[i, iter]
            self.Edpp[i, iter] = self.Vd[i, iter] + gen.Ra*self.Id[i, iter] + gen.Xqpp*self.Iq[i, iter]

    def __state_equations(self, x, y, iter = 0):
        """Calculates the state equations of a synchronous generator.

        Args:
            iter: Current iteration index

        Returns:
            A NumPy array containing the derivatives of omega, delta, and the fluxes, for each generator

        """

        d_omega_dt = np.array([
            ((0.5/gen.H) * (self.Pmgap[i,iter] - self.Pegap[i,iter]) for i, gen in enumerate(self.system.generators))
        ])
        d_delta_dt = self.w[:,iter]*self.system.we
        d_Edp_dt = np.array([
            ((1/gen.Tq0p) * (-self.Edp[i,iter] - (gen.Xq - gen.Xqp)*self.Iq[i, iter]) for i, gen in enumerate(self.system.generators))
        ])
        d_Eqp_dt = np.array([
            ((1/gen.Td0p) * (-self.Eqp[i,iter] - (gen.Xd - gen.Xdp)*self.Id[i, iter] + self.Vexc[i,iter]) for i, gen in enumerate(self.system.generators))
        ])
        d_Edpp_dt = np.array([
            ((1/gen.Tq0pp) * (-self.Edpp[i,iter] - (gen.Xqp - gen.Xqpp)*self.Iq[i, iter] + self.Edp[i,iter]) for i, gen in enumerate(self.system.generators))
        ])
        d_Eqpp_dt = np.array([
            ((1/gen.Td0pp) * (-self.Eqpp[i,iter] - (gen.Xdp - gen.Xdpp)*self.Id[i, iter] + self.Eqp[i,iter]) for i, gen in enumerate(self.system.generators))
        ])

        return np.array([d_omega_dt, d_delta_dt, d_Edp_dt, d_Eqp_dt, d_Edpp_dt, d_Eqpp_dt]).reshape((6, self.system.n_gens))

    def solve(self, method = 'RK4'):
        # Compute initial values
        self.__initialize_values()

        # Get explicit solver
        solver = RK4

        # Now, start simulation
        for i in range(self.N-1):

            # Compute normal equations
            self.__compute_non_dynamic_values(i)

            y_old = self.__state_equations(self.tspan[i], y_old, self.__dt, i)
            y_new = solver.step(self.__state_equations, self.tspan[i], y_old, self.__dt, i)

            # Assin new values
            self.omega[:,i+1] = y_new[0,:]
            self.delta[:,i+1] = y_new[1,:]
            self.Edp[:,i+1] = y_new[2,:]
            self.Eqp[:,i+1] = y_new[3, :]
            self.Edpp[:, i+1] = y_new[4, :]
            self.Eqpp[:, i+1] = y_new[5, :]

    def __create_state_equations_model(self):
        m = GEKKO()

        pass

    def __create_state_equations_model_variables(self, m):
        self.wn = m.Array(m.Var, dim = (self.system.n_gens,))
        self.dn = m.Array(m.Var, dim = (self.system.n_gens,))

    def __solve_state_equations(self):
        pass


    @property
    def N(self):
        return len(self.tspan)





