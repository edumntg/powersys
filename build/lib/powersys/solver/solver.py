from ..models.powersystem import PowerSystem
import numpy as np
import pandas as pd
from gekko import GEKKO

class Solver:
    def __init__(self, model: PowerSystem):
        self.model = model
        self.solved = False

    def solve(self, display = True):
        pass

    def extract_results(self):
        if not self.solved:
            raise "Model not solved yet. Call solve()"
        
        # Create dataframes for results
        data = []
        for bus in self.model.buses:
            # Append results as: V, theta, Pload, Qload, Pgen, Qgen
            row = [bus.id, self.V[bus.id][0], self.theta[bus.id][0], bus.Pload, bus.Qload]
            Pgen = self.Pgen_fixed[bus.id][0] + np.sum([self.Pgen[gen.bus][0] for gen in self.model.generators if gen.bus == bus.id])
            Qgen = self.Qgen_fixed[bus.id][0] + np.sum([self.Qgen[gen.bus][0] for gen in self.model.generators if gen.bus == bus.id])
            row += [Pgen, Qgen, bus.Vmin, bus.Vmax]

            data.append(row)
        
        # Create DataFrame
        bus_df = pd.DataFrame(data, columns = ['id', 'V', 'theta', 'Pload', 'Qload', 'Pgen', 'Qgen', 'Vmin', 'Vmax'])

        # Create generators DataFrame
        data = []
        for gen in self.model.generators:
            row = [gen.id, self.lg[gen.id][0], self.Pgen[gen.id][0], self.Qgen[gen.id][0], gen.cost(self.Pgen[gen.id][0]), gen.Pmin, gen.Pmax, gen.Qmin, gen.Qmax]
            data.append(row)
        
        gen_df = pd.DataFrame(data, columns = ['id', 'active', 'Pgen', 'Qgen', 'cost', 'Pmin', 'Pmax', 'Qmin', 'Qmax'])

        #Create lines DataFrame
        data = []
        for line in self.model.lines:
            i = line.from_bus
            k = line.to_bus
            row = [line.id, i, k, self.l[line.id].value[0], self.Pflow[i,k][0], self.Pflow[k,i][0], self.Qflow[i,k][0], self.Qflow[k,i][0], line.S(self.Pflow[i,k][0], self.Qflow[i,k][0]), line.Plosses(self.Pflow[i,k][0], self.Pflow[k,i][0]), line.mva]
            data.append(row)
        
        lines_df = pd.DataFrame(data, columns = ['id', 'from', 'to', 'active', 'Pout', 'Pin', 'Qout', 'Qin', 'S', 'Ploss', 'mva'])

        return bus_df, gen_df, lines_df
    
    def construct_model_solver(self):
        solver = GEKKO()

        # Voltage variables
        self.V = solver.Array(solver.Var, dim = (self.model.n_buses,), value = 1.0)

        # Angle variables
        self.theta = solver.Array(solver.Var, dim = (self.model.n_buses,), lb = -np.pi, ub = np.pi, value = 0.0)

        # Active power generated
        self.Pgen = solver.Array(solver.Var, dim = (self.model.n_gens,))
        self.Pgen_fixed = solver.Array(solver.Param, dim = (self.model.n_buses,))
        for bus in self.model.buses:
            self.Pgen_fixed[bus.id].VALUE = bus.Pgen_fixed

        # Reactive power generated
        self.Qgen = solver.Array(solver.Var, dim = (self.model.n_gens,))
        self.Qgen_fixed = solver.Array(solver.Param, dim = (self.model.n_buses,))
        for bus in self.model.buses:
            self.Qgen_fixed[bus.id].VALUE = bus.Qgen_fixed
        
        # Flow through lines
        self.Pflow = solver.Array(solver.Var, dim = (self.model.n_buses, self.model.n_buses))
        self.Qflow = solver.Array(solver.Var, dim = (self.model.n_buses, self.model.n_buses))

        # On/Off for lines
        self.l = solver.Array(solver.Var, dim = (self.model.n_lines,), lb = 0, ub = 1, value = 1, integer = True)

        # On/Off for generators
        self.lg = solver.Array(solver.Var, dim = (self.model.n_gens,), lb = 0, ub = 1, value = 1, integer = True)

        self.construct_optim_constr(solver)

        return solver
    
    def construct_optim_constr(self, m):
        # Kirchoff law for active power
        m.Equations([
            self.optim_constr_kirchoff_P(bus) for bus in self.model.buses
        ])
        
        # Kirchoff law for reactive power
        m.Equations([
            self.optim_constr_kirchoff_Q(bus) for bus in self.model.buses
        ])

        # Equations for active-power flow through lines
        m.Equations([
            self.optim_constr_line_P_fromto(line, m) for line in self.model.lines
        ])
        m.Equations([
            self.optim_constr_line_P_tofrom(line, m) for line in self.model.lines
        ])

        # # Equations for reactive-power flow through lines
        m.Equations([
            self.optim_constr_line_Q_fromto(line, m) for line in self.model.lines
        ])
        m.Equations([
            self.optim_constr_line_Q_tofrom(line, m) for line in self.model.lines
        ])

        # Add constraints for reference buses
        m.Equations([
            self.optim_constr_reference_bus_angle(bus) for bus in self.model.buses if bus.reference
        ])
    
    def optim_constr_kirchoff_P(self, bus):
        Pflow = 0.0
        Pgen = self.Pgen_fixed[bus.id]

        # Find if there is a generator connected at this bus
        for gen in self.model.generators:
            if gen.bus == bus.id:
                Pgen += self.Pgen[gen.id]

        for line in self.model.lines:
            if line.from_bus == bus.id:
                Pflow += self.Pflow[line.from_bus, line.to_bus]

        for line in self.model.lines:
            if line.to_bus == bus.id:
                Pflow += self.Pflow[line.to_bus, line.from_bus]

        return Pgen == bus.Pload + Pflow
    
    def optim_constr_kirchoff_Q(self, bus):
        Qflow = 0.0
        Qgen = self.Qgen_fixed[bus.id]

        # Find if there is a generator connected at this bus
        for gen in self.model.generators:
            if gen.bus == bus.id:
                Qgen += self.Qgen[gen.id]

        for line in self.model.lines:
            if line.from_bus == bus.id:
                Qflow += self.Qflow[line.from_bus, line.to_bus]

        for line in self.model.lines:
            if line.to_bus == bus.id:
                Qflow += self.Qflow[line.to_bus, line.from_bus]

        return Qgen == bus.Qload + Qflow
    
    def optim_constr_line_P_fromto(self, line, m):
        i = line.from_bus
        k = line.to_bus

        return self.Pflow[i,k] == self.l[line.id]*((-self.model.G[i,k] + self.model.g[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(self.model.G[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.model.B[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def optim_constr_line_P_tofrom(self, line, m):
        i = line.to_bus
        k = line.from_bus

        return self.Pflow[i,k] == self.l[line.id]*((-self.model.G[i,k] + self.model.g[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(self.model.G[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.model.B[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def optim_constr_line_Q_fromto(self, line, m):
        i = line.from_bus
        k = line.to_bus

        return self.Qflow[i,k] == self.l[line.id]*((self.model.B[i,k] - self.model.b[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(-self.model.B[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.model.G[i,k]*m.sin(self.theta[i] - self.theta[k])))
    
    def optim_constr_line_Q_tofrom(self, line, m):
        i = line.to_bus
        k = line.from_bus

        return self.Qflow[i,k] == self.l[line.id]*((self.model.B[i,k] - self.model.b[i,k])*self.V[i]**2 + self.V[i]*self.V[k]*(-self.model.B[i,k]*m.cos(self.theta[i] - self.theta[k]) + self.model.G[i,k]*m.sin(self.theta[i] - self.theta[k])))

    def optim_constr_reference_bus_angle(self, bus):
        return self.theta[bus.id] == 0.0
    
    def loadflow_constr_reference_bus_voltage(self, bus):
        return self.V[bus.id] == bus.V

    def assign_results(self):
        # Assign results
        for bus in self.model.buses:
            bus.V = self.V[bus.id][0]
            bus.theta = self.theta[bus.id][0]
            bus.Pgen = np.sum([self.Pgen[gen.id][0] for gen in self.model.generators if gen.bus == bus.id])
            bus.Qgen = np.sum([self.Qgen[gen.id][0] for gen in self.model.generators if gen.bus == bus.id])

        for gen in self.model.generators:
            gen.Pgen = self.Pgen[gen.id][0]
            gen.Qgen = self.Qgen[gen.id][0]