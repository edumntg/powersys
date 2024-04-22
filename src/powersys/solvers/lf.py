from ..models.powersystem import PowerSystem
from .solver import *

class LF(Solver):

    def __init__(self, model: PowerSystem):
        super().__init__(model)

    def solve(self, disp = True, method = "default"):
        print("METHOD:", method)
        if not self.model:
            raise "No PowerSystem object declared"
        
        # Construct Ybus
        if self.model.Ybus is None:
            self.model.construct_ybus()

        # Construct Optim model
        if method == "default":
            solver_model = self.construct_model_solver()
        elif method == "gauss-seidel" or method == "gs":
            solver_model = self.construct_iterative_solver(self.model, "gauss-seidel")
        elif method =="newton-raphson" or method == "nr":
            #solver_model = self.construct_newton_raphson_model()
            pass
        elif method == "scipy" or method == "fsolve":
            solver_model = self.construct_iterative_solver(self.model, "scipy")
        else:
            raise "Invalid LF solver method specified. Got: " + method

        # Now, solve
        if method == "default":
            solver_model.solve(disp = disp)

            V = np.array([self.V[bus.id][0] for bus in self.model.buses])
            theta = np.array([self.theta[bus.id][0] for bus in self.model.buses])
            Pgen = np.array([self.Pgen[gen.id][0] for gen in self.model.generators])
            Qgen = np.array([self.Qgen[gen.id][0] for gen in self.model.generators])
        else:
            V, theta, Pgen, Qgen = solver_model.solve(disp = disp)
        
        self.solved = True
        self.model.lf_solved = True

        # Assign results
        self.assign_results(V, theta, Pgen, Qgen)

    def construct_model_solver(self):
        solver = super().construct_model_solver()

        self.__construct_load_flow_constraints(solver)
        
        return solver

    def __construct_load_flow_constraints(self, m):

        # Set voltage magnitude fixed for PV buses
        # NOTE: For slack bus, the constraint for volt magnitude has already been set in the Solver instance. Same for angle
        m.Equations([
            self.optim_constr_fixed_bus_magnitude(bus) for bus in self.model.buses if bus.type == 2
        ])

        # Set lines always active
        m.Equations([
            self.l[line.id] == 1 for line in self.model.lines
        ])

        # Set generators always active
        m.Equations([
            self.lg[gen.id] == 1 for gen in self.model.generators
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
            if not self.model.get_bus(gen.bus).type == 1:
                m.Equation(self.Pgen[gen.id] == 0.0)