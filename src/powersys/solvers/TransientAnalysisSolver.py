from ..models.PowerSystem import PowerSystem
from ..math.explicit.RK4 import RK4
import numpy as np
from gekko import GEKKO

class TransientAnalysisSolver(object):

    def __init__(self, system: PowerSystem, tspan: np.array):
        self.system = system
        self.solved = False
        self.tspan = tspan

        self.__initialize_values() # Initialize all vectors and pre-event values

    def __initialize_values(self):
        if not self.system:
            raise "No PowerSystem object given"
        
        N = len(self.tspan) # Number of time points

        self.d = np.zeros((self.system.n_gens, N))
        self.w = np.zeros((self.system.n_gens, N))
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

        self.ev = np.zeros((self.system.n_gens, 1))

        # Assign initial values
        self.d[:,0] = np.array([gen.df for gen in self.system.generators])
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

    def __state_equations(self):
        # For each generator, we will have 7 equations
        for i, gen in enumerate(self.system.generators):
            self.dw[i] = (0.5/gen.H) * (self.Pmgap[i] - self.Pegap[i])
            self.dd[i] = self.w[i]*self.system.we
            self.dVmed[i] = (1.0/gen.Tmed)*(gen.Kmed*np.abs(self.Vt[i] - self.Vmed[i]))
            if gen.USE_PSS:
                self.ev[i] = self.V[i] - self.Vmed[i] + self.Vpc2[i]
            else:
                self.ev[i] = self.V[i] - self.Vmed[i]

            self.dVvi[i] = gen.Kv*self.ev[i]
            
            # PSS
            self.dVw[i] = gen.Kest*self.dw[i] - self.Vw[i]/gen.Tw
            self.dVpc1[i] = (1/gen.T2)*(self.Vw[i] + gen.T1*self.dVw[i] - self.Vpc1[i])
            self.dVpc2[i] = (1/gen.T4)*(self.Vpc1[i] + gen.T3*self.dVpc1[i] - self.Vpc2[i])

            self.Vvp[i] = gen.Kp*self.ev[i]

            if gen.USE_PSS:
                self.Vvd[i] = gen.Kd*(-self.dVmed[i] + self.dVpc2[i])
            else:
                self.Vvd[i] = -gen.Kd*self.dVmed[i]

            self.Vx[i] = self.Vvp[i] + self.Vvd[i] + self.Vvi[i]
            self.dVa[i] = (1.0/gen.Ta)*(gen.Ka*self.Vx[i] - self.Va[i])
            self.dVexc[i] = (1.0/gen.Texc)*(gen.Kexc*self.Va[i] - self.Vexc[i])

            self.dEqp[i] = (1.0/gen.Td0p)*(-self.Eqp[i] + (gen.Xd - gen.Xdp)*self.Id[i] + self.Vexc[i])
            self.dEdp[i] = (1.0/gen.Tq0p)*(-self.Edp[i] - (gen.Xq - gen.Xqp)*self.Iq[i])
            self.dEqpp[i] = (1.0/gen.Td0pp)*(-self.Eqpp[i] + (gen.Xd - gen.Xdpp)*self.Id[i] + self.Eqp[i])
            self.dEdpp[i] = (1.0/gen.Tq0pp)*(-self.Edpp[i] - (gen.Xq - gen.Xqpp)*self.Iq[i] + self.Edp[i])

            self.Xp[i] = -self.w[i]/gen.R

            self.dXv[i] = (1.0/gen.Tv)*(gen.Kv*(self.Pc[i] + gen.Xp) - self.Xv[i])
            self.dPc[i] = -gen.Ki*self.w[i]


    def solve(self, method = 'RK4'):
        # Compute initial values
        self.__initialize_values()

        # Get explicit solver
        solver = RK4

        # Now, start simulation
        for i in range(self.N-1):

            pass

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





