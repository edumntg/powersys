from ..models.PowerSystem import PowerSystem
from gekko import GEKKO
import pandas as pd
import numpy as np

class OPFSolver(object):

    def __init__(self, system: PowerSystem):
        self.system = system
        self.opf_solved = False
        self.opf_solution = None
        self.opf_solver = None

    def solve(self, disp = True):
        if not self.system:
            raise "No PowerSystem object declared"
        
        # Construct Ybus
        self.system.construct_ybus()

        # Construct Pyomo model
        model = self.__construct_optim_model()

        # Now, solve
        model.solve(disp = disp)
        self.opf_solved = True

        # Assign results
        for bus in self.system.buses:
            bus.V = self.V[bus.id][0]
            bus.theta = self.theta[bus.id][0]
            bus.Pgen = self.Pgen_fixed[bus.id][0] + np.sum([self.Pgen[gen.id][0] for gen in self.system.generators if gen.bus == bus.id])
            bus.Qgen = self.Qgen_fixed[bus.id][0] + np.sum([self.Qgen[gen.id][0] for gen in self.system.generators if gen.bus == bus.id])

        for gen in self.system.generators:
            gen.Pgen = self.Pgen[gen.id][0]
            gen.Qgen = self.Qgen[gen.id][0]

        # Load lagrange and kuhn-tucker multipliers
        # multipliers = np.loadtxt(model.path + '/apm_lam.txt')
        # self.multipliers = multipliers

    def extract_results(self):
        if not self.opf_solved and not self.load_flow_solved:
            raise "OPF not solved yet"
        
        # Create dataframes for results
        data = []
        for bus in self.system.buses:
            # Append results as: V, theta, Pload, Qload, Pgen, Qgen
            row = [bus.id, self.V[bus.id][0], self.theta[bus.id][0], bus.Pload, bus.Qload]
            Pgen = self.Pgen_fixed[bus.id][0] + np.sum([self.Pgen[gen.bus][0] for gen in self.system.generators if gen.bus == bus.id])
            Qgen = self.Qgen_fixed[bus.id][0] + np.sum([self.Qgen[gen.bus][0] for gen in self.system.generators if gen.bus == bus.id])
            row += [Pgen, Qgen, bus.Vmin, bus.Vmax]

            data.append(row)
        
        # Create DataFrame
        bus_df = pd.DataFrame(data, columns = ['id', 'V', 'theta', 'Pload', 'Qload', 'Pgen', 'Qgen', 'Vmin', 'Vmax'])

        # Create generators DataFrame
        data = []
        for gen in self.system.generators:
            row = [gen.id, self.lg[gen.id][0], self.Pgen[gen.id][0], self.Qgen[gen.id][0], gen.cost(self.Pgen[gen.id][0]), gen.Pmin, gen.Pmax, gen.Qmin, gen.Qmax]
            data.append(row)
        
        gen_df = pd.DataFrame(data, columns = ['id', 'active', 'Pgen', 'Qgen', 'cost', 'Pmin', 'Pmax', 'Qmin', 'Qmax'])

        #Create lines DataFrame
        data = []
        for line in self.system.lines:
            i = line.from_bus
            k = line.to_bus
            row = [line.id, i, k, self.l[line.id].value[0], self.Pflow[i,k][0], self.Pflow[k,i][0], self.Qflow[i,k][0], self.Qflow[k,i][0], line.S(self.Pflow[i,k][0], self.Qflow[i,k][0]), line.Plosses(self.Pflow[i,k][0], self.Pflow[k,i][0]), line.mva]
            data.append(row)
        
        lines_df = pd.DataFrame(data, columns = ['id', 'from', 'to', 'active', 'Pout', 'Pin', 'Qout', 'Qin', 'S', 'Ploss', 'mva'])

        return bus_df, gen_df, lines_df
    
    def __opf_objective(self):
        total_cost = np.sum([gen.cost(self.Pgen[gen.id]) for gen in self.system.generators])

        return total_cost

    def __construct_optim_constr(self, m):
        # Kirchoff law for active power
        m.Equations([
            self.__opf_constr_kirchoff_P(bus) for bus in self.system.buses
        ])
        
        # Kirchoff law for reactive power
        m.Equations([
            self.__opf_constr_kirchoff_Q(bus) for bus in self.system.buses
        ])

        # Equations for active-power flow through lines
        m.Equations([
            self.__opf_constr_line_P_fromto(line, m) for line in self.system.lines
        ])
        m.Equations([
            self.__opf_constr_line_P_tofrom(line, m) for line in self.system.lines
        ])

        # # Equations for reactive-power flow through lines
        m.Equations([
            self.__opf_constr_line_Q_fromto(line, m) for line in self.system.lines
        ])
        m.Equations([
            self.__opf_constr_line_Q_tofrom(line, m) for line in self.system.lines
        ])

        # Add constraints for reference buses
        m.Equations([
            self.__opf_constr_reference_bus_angle(bus) for bus in self.system.buses if bus.reference
        ])

    def __construct_opf_contraints(self, m):
        # Bus voltage limits
        m.Equations([
            self.__opf_constr_bus_voltage_min(bus) for bus in self.system.buses
        ])
        m.Equations([
            self.__opf_constr_bus_voltage_max(bus) for bus in self.system.buses
        ])

        # Min active power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Pmin(gen) for gen in self.system.generators
        ])

        # Max active power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Pmax(gen) for gen in self.system.generators
        ])

        # Min reactive power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Qmin(gen) for gen in self.system.generators
        ])

        # Max reactive power generated by each generator
        m.Equations([
            self.__opf_constr_gen_Qmax(gen) for gen in self.system.generators
        ])

        # Transmission lines limits
        m.Equations([
            self.__opf_constr_line_max_mva_fromto(line) for line in self.system.lines
        ])
        m.Equations([
            self.__opf_constr_line_max_mva_tofrom(line) for line in self.system.lines
        ])
    
    def __construct_optim_model(self):
        m = GEKKO()

        m.options.SOLVER = 1 # apopt for MINLP
            # m.options.DIAGLEVEL = 2 # for multipliers
        # Voltage variables
        self.V = m.Array(m.Var, dim = (self.system.n_buses,), value = 1.0)

        # Angle variables
        self.theta = m.Array(m.Var, dim = (self.system.n_buses,), lb = -np.pi, ub = np.pi, value = 0.0)

        # Active power generated
        self.Pgen = m.Array(m.Var, dim = (self.system.n_gens,))
        self.Pgen_fixed = m.Array(m.Param, dim = (self.system.n_buses,))
        for bus in self.system.buses:
            self.Pgen_fixed[bus.id].VALUE = bus.Pgen_fixed

        # Reactive power generated
        self.Qgen = m.Array(m.Var, dim = (self.system.n_gens,))
        self.Qgen_fixed = m.Array(m.Param, dim = (self.system.n_buses,))
        for bus in self.system.buses:
            self.Qgen_fixed[bus.id].VALUE = bus.Qgen_fixed
        
        # Flow through lines
        self.Pflow = m.Array(m.Var, dim = (self.system.n_buses, self.system.n_buses))
        self.Qflow = m.Array(m.Var, dim = (self.system.n_buses, self.system.n_buses))

        # On/Off for lines
        self.l = m.Array(m.Var, dim = (self.system.n_lines,), lb = 0, ub = 1, value = 1, integer = True)

        # On/Off for generators
        self.lg = m.Array(m.Var, dim = (self.system.n_gens,), lb = 0, ub = 1, value = 1, integer = True)

        # Create objective function
        m.Obj(self.__opf_objective())

        # Construct constraints which are common for both opf/loadflow
        self.__construct_optim_constr(m)

        self.__construct_opf_contraints(m)

        self.m = m

        return self.m

    def __opf_constr_gen_Pmin(self, gen):
        return self.Pgen[gen.id] >= self.lg[gen.id]*gen.Pmin
    
    def __opf_constr_gen_Pmax(self, gen):
        return self.Pgen[gen.id] <= self.lg[gen.id]*gen.Pmax
    
    def __opf_constr_gen_Qmin(self, gen):
        return self.Qgen[gen.id] >= self.lg[gen.id]*gen.Qmin
    
    def __opf_constr_gen_Qmax(self, gen):
        return self.Qgen[gen.id] <= self.lg[gen.id]*gen.Qmax
    
    def __opf_constr_kirchoff_P(self, bus):
        Pflow = 0.0
        Pgen = self.Pgen_fixed[bus.id]

        # Find if there is a generator connected at this bus
        for gen in self.system.generators:
            if gen.bus == bus.id:
                Pgen += self.Pgen[gen.id]

        for line in self.system.lines:
            if line.from_bus == bus.id:
                Pflow += self.Pflow[line.from_bus, line.to_bus]

        for line in self.system.lines:
            if line.to_bus == bus.id:
                Pflow += self.Pflow[line.to_bus, line.from_bus]

        return Pgen == bus.Pload + Pflow
    
    def __opf_constr_kirchoff_Q(self, bus):
        Qflow = 0.0
        Qgen = self.Qgen_fixed[bus.id]

        # Find if there is a generator connected at this bus
        for gen in self.system.generators:
            if gen.bus == bus.id:
                Qgen += self.Qgen[gen.id]

        for line in self.system.lines:
            if line.from_bus == bus.id:
                Qflow += self.Qflow[line.from_bus, line.to_bus]

        for line in self.system.lines:
            if line.to_bus == bus.id:
                Qflow += self.Qflow[line.to_bus, line.from_bus]

        return Qgen == bus.Qload + Qflow
    
    def __opf_constr_line_P_fromto(self, line, m):
        i = line.from_bus
        k = line.to_bus

        return self.Pflow[i,k] == self.l[line.id]*((-self.system.G[i,k] + self.system.g[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(self.system.G[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.system.B[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def __opf_constr_line_P_tofrom(self, line, m):
        i = line.to_bus
        k = line.from_bus

        return self.Pflow[i,k] == self.l[line.id]*((-self.system.G[i,k] + self.system.g[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(self.system.G[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.system.B[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def __opf_constr_line_Q_fromto(self, line, m):
        i = line.from_bus
        k = line.to_bus

        return self.Qflow[i,k] == self.l[line.id]*((self.system.B[i,k] - self.system.b[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(-self.system.B[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.system.G[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def __opf_constr_line_Q_tofrom(self, line, m):
        i = line.to_bus
        k = line.from_bus

        return self.Qflow[i,k] == self.l[line.id]*((self.system.B[i,k] - self.system.b[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(-self.system.B[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.system.G[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def __opf_constr_line_max_mva_fromto(self, line):
        # Get apparent power

        #return self.Pflow[line.from_bus, line.to_bus]**2 + self.Qflow[line.from_bus, line.to_bus]**2 <= line.mva**2
        return self.Pflow[line.from_bus, line.to_bus]**2 <= line.mva**2
    
    def __opf_constr_line_max_mva_tofrom(self, line):
        # Get apparent power

        #return self.Pflow[line.to_bus, line.from_bus]**2 + self.Qflow[line.to_bus, line.from_bus]**2 <= line.mva**2
        return self.Pflow[line.to_bus, line.from_bus]**2 <= line.mva**2

    def __opf_constr_reference_bus_angle(self, bus):
        return self.theta[bus.id] == 0.0
    
    def __opf_constr_bus_voltage_min(self, bus):
        return self.V[bus.id] >= bus.Vmin
    
    def __opf_constr_bus_voltage_max(self, bus):
        return self.V[bus.id] <= bus.Vmax
        