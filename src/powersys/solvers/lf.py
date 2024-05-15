from ..models.powersystem import PowerSystem
from .solver import *
from ..math.iterative import IterativeArgs

class LF(Solver):

    def __init__(self, model: PowerSystem):
        super().__init__(model)

    def solve(self, disp = True, method = "default", **kwargs):
        print("METHOD:", method)
        if not self.model:
            raise "No PowerSystem object declared"
        
        # Construct Ybus
        if self.model.Ybus is None:
            self.model.construct_ybus()

        # Construct Optim model
        if method == "default":
            self.construct_model_solver()
        elif method == "gauss-seidel" or method == "gs":
            self.construct_iterative_solver(self.model, "gauss-seidel", IterativeArgs(**kwargs))
        elif method =="newton-raphson" or method == "nr":
            self.construct_iterative_solver(self.model, "newton-raphson", IterativeArgs(**kwargs))
        elif method == "scipy" or method == "fsolve":
            self.construct_iterative_solver(self.model, "scipy", IterativeArgs(**kwargs))
        else:
            raise "Invalid LF solver method specified. Got: " + method

        # Now, solve
        solver_model, variables = self.state_dict()['solver'], self.state_dict()['variables']
        if method == "default":
            solver_model.solve(disp = disp)
            V = np.array([variables['V'][bus.id][0] for bus in self.model.buses])
            theta = np.array([variables['theta'][bus.id][0] for bus in self.model.buses])
            Pgen = np.array([variables['Pgen'][gen.id][0] for gen in self.model.generators])
            Qgen = np.array([variables['Qgen'][gen.id][0] for gen in self.model.generators])

            # Assign results to system
            for bus in self.model.buses:
                bus.V = V[bus.id]
                bus.theta = theta[bus.id]
            
            for gen in self.model.generators:
                gen.Pgen = Pgen[gen.id]
                gen.Qgen = Qgen[gen.id]
        else:
            V, theta, Pgen, Qgen = solver_model.solve(self.state_dict(), disp = disp) # here, Pgen is for each bus
            
            # Assign results
            for bus in self.model.buses:
                bus.V = V[bus.id]
                bus.theta = theta[bus.id]

                gen = self.model.get_gen_by_bus(bus.id)
                if gen:
                    gen.Pgen = Pgen[bus.id] - bus.Pgen_fixed
                    gen.Qgen = Qgen[bus.id] - bus.Qgen_fixed

            Pgenbus = np.array([[Pgen[bus.id] - bus.Pgen_fixed] for bus in self.model.buses if self.model.get_gen_by_bus(bus.id)])
            Qgenbus = np.array([[Qgen[bus.id] - bus.Qgen_fixed] for bus in self.model.buses if self.model.get_gen_by_bus(bus.id)])
            variables = {
                'V': V,
                'theta': theta,
                'Pgen': Pgenbus,
                'Qgen': Qgenbus,
                'lg': np.array([[gen.active] for gen in self.model.generators]),
                'l': np.array([[1] for line in self.model.lines])
            }

            # Compute flows
            Pflow = np.zeros((self.model.N, self.model.N, 1))
            Qflow = np.zeros((self.model.N, self.model.N, 1))
            for line in self.model.lines:
                i = line.from_bus
                k = line.to_bus
                Pflow[i,k,0] = ((-self.model.G[i,k] + self.model.g[i,k])*V[i]**2 + V[i]*V[k]*(self.model.G[i,k]*np.cos(theta[i] - theta[k]) + self.model.B[i,k]*np.sin(theta[i] - theta[k])))
                Pflow[k,i,0] = ((-self.model.G[k,i] + self.model.g[k,i])*V[k]**2 + V[k]*V[i]*(self.model.G[k,i]*np.cos(theta[k] - theta[i]) + self.model.B[k,i]*np.sin(theta[k] - theta[i])))

                Qflow[i,k,0] = ((self.model.B[i,k] - self.model.b[i,k])*V[i]**2 + V[i]*V[k]*(-self.model.B[i,k]*np.cos(theta[i] - theta[k]) + self.model.G[i,k]*np.sin(theta[i] - theta[k])))
                Qflow[k,i,0] = ((self.model.B[k,i] - self.model.b[k,i])*V[k]**2 + V[k]*V[i]*(-self.model.B[k,i]*np.cos(theta[k] - theta[i]) + self.model.G[k,i]*np.sin(theta[k] - theta[i])))
            
            variables['Pflow'] = Pflow
            variables['Qflow'] = Qflow
            self.set_variables(variables)

        self.solved = True
        self.model.lf_solved = True

        # Assign results
        self.assign_results(V, theta, Pgen, Qgen)

    def construct_model_solver(self):
        super().construct_model_solver()

        self.__construct_load_flow_constraints()

    def __construct_load_flow_constraints(self):

        m, variables = self.state_dict()['solver'], self.state_dict()['variables']
        # Set voltage magnitude fixed for PV buses
        # NOTE: For slack bus, the constraint for volt magnitude has already been set in the Solver instance. Same for angle
        m.Equations([
            self.optim_constr_fixed_bus_magnitude(bus) for bus in self.model.buses if bus.type == PowerSystem.PV
        ])

        # Set lines always active
        m.Equations([
            variables['l'][line.id] == 1 for line in self.model.lines
        ])

        # Set generators always active
        m.Equations([
            variables['lg'][gen.id] == 1 for gen in self.model.generators
        ])

        # For all buses containing a generator, we assume that bus is PV
        # for bus in self.model.buses:
        #     is_pv = False
        #     for gen in self.model.generators:
        #         if gen.bus == bus.id:
        #             is_pv = True
        #             break
            
        #     if is_pv and not bus.type == 1:
        #         m.Equation(self.V[bus.id] == bus.V)

        # Now, for generators, set P = 0 (so they only supply Q). This is because PV buses have already a P declared that does not changes
        for gen in self.model.generators:
            if not self.model.get_bus(gen.bus).type == PowerSystem.SLACK:
                m.Equation(variables['Pgen'][gen.id] == 0.0)