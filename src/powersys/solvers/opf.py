from ..models.powersystem import PowerSystem
import numpy as np
from ..solvers.solver import Solver 

class OPF(Solver):

    def __init__(self, model: PowerSystem):
        super().__init__(model)

    def solve(self, disp = True):
        if not self.model:
            raise "No PowerSystem object declared"
        
        # Construct Ybus
        self.model.construct_ybus()

        # Construct Pyomo model
        solver_model = self.construct_model_solver()

        # Now, solve
        solver_model.solve(disp = disp)
        self.solved = True

        # Assign results
        self.assign_results()

        # Load lagrange and kuhn-tucker multipliers
        # multipliers = np.loadtxt(model.path + '/apm_lam.txt')
        # self.multipliers = multipliers

    def construct_model_solver(self):

        solver = super().construct_model_solver()

        # Create objective function
        solver.Obj(self.__opf_objective())

        self.__construct_opf_contraints(solver)

        return solver
    
    def __opf_objective(self):
        total_cost = np.sum([gen.cost(self.Pgen[gen.id]) for gen in self.model.generators])

        return total_cost

    def __construct_opf_contraints(self, m):
        # Bus voltage limits
        m.Equations([
            self.__opf_constr_bus_voltage_min(bus) for bus in self.model.buses
        ])
        m.Equations([
            self.__opf_constr_bus_voltage_max(bus) for bus in self.model.buses
        ])

        # Min active power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Pmin(gen) for gen in self.model.generators
        ])

        # Max active power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Pmax(gen) for gen in self.model.generators
        ])

        # Min reactive power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Qmin(gen) for gen in self.model.generators
        ])

        # Max reactive power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Qmax(gen) for gen in self.model.generators
        ])

        # Transmission lines limits
        m.Equations([
            self.__opf_constr_line_max_mva_fromto(line) for line in self.model.lines
        ])
        m.Equations([
            self.__opf_constr_line_max_mva_tofrom(line) for line in self.model.lines
        ])
    
    def __opf_constr_gen_Pmin(self, gen):
        return self.Pgen[gen.id] >= self.lg[gen.id]*gen.Pmin
    
    def __opf_constr_gen_Pmax(self, gen):
        return self.Pgen[gen.id] <= self.lg[gen.id]*gen.Pmax
    
    def __opf_constr_gen_Qmin(self, gen):
        return self.Qgen[gen.id] >= self.lg[gen.id]*gen.Qmin
    
    def __opf_constr_gen_Qmax(self, gen):
        return self.Qgen[gen.id] <= self.lg[gen.id]*gen.Qmax
    
        i = line.to_bus
        k = line.from_bus

        return self.Qflow[i,k] == self.l[line.id]*((self.model.B[i,k] - self.model.b[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(-self.model.B[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.model.G[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def __opf_constr_line_max_mva_fromto(self, line):
        # Get apparent power

        #return self.Pflow[line.from_bus, line.to_bus]**2 + self.Qflow[line.from_bus, line.to_bus]**2 <= line.mva**2
        return self.Pflow[line.from_bus, line.to_bus]**2 <= line.mva**2
    
    def __opf_constr_line_max_mva_tofrom(self, line):
        # Get apparent power

        #return self.Pflow[line.to_bus, line.from_bus]**2 + self.Qflow[line.to_bus, line.from_bus]**2 <= line.mva**2
        return self.Pflow[line.to_bus, line.from_bus]**2 <= line.mva**2

    def __opf_constr_bus_voltage_min(self, bus):
        return self.V[bus.id] >= bus.Vmin
    
    def __opf_constr_bus_voltage_max(self, bus):
        return self.V[bus.id] <= bus.Vmax
        