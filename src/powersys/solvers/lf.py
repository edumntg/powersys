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
            solver_model = self.construct_gauss_seidel_model(self.model)
        elif method =="newton-raphson" or method == "nr":
            #solver_model = self.construct_newton_raphson_model()
            pass
        else:
            raise "Invalid LF solver method specified. Got: " + method
        
        # SCIPY TEST
        func, x0 = self.construct_scipy_fsolve(self.model)
        print(self.solve_scipy_fsolve(func, x0, self.model))

        # Now, solve
        solver_model.solve(disp = disp)
        self.solved = True
        self.model.lf_solved = True
        
        # Assign results
        self.assign_results()

    def construct_model_solver(self):
        solver = super().construct_model_solver()

        self.__construct_load_flow_constraints(solver)
        
        return solver

    def __construct_load_flow_constraints(self, m):
        m.Equations([
            self.__loadflow_constr_reference_bus_voltage(bus) for bus in self.model.buses if bus.reference
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
        for bus in self.model.buses:
            is_pv = False
            for gen in self.model.generators:
                if gen.bus == bus.id:
                    is_pv = True
                    break
            
            if is_pv and not bus.reference:
                m.Equation(self.V[bus.id] == bus.V)

        # Now, for generators, set P = 0 (so they only supply Q)
        for gen in self.model.generators:
            if not self.model.get_bus(gen.bus).reference:
                m.Equation(self.Pgen[gen.id] == 0.0)

    def __loadflow_constr_reference_bus_voltage(self, bus):
        return self.V[bus.id] == bus.V

    def construct_gauss_seidel_model(self, model: PowerSystem):
        solver = super().construct_iterative_solver(model, "gauss-seidel")

        return solver