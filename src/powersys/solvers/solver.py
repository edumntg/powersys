from ..models.powersystem import PowerSystem
from ..math.iterative import *
import numpy as np
import pandas as pd
from gekko import GEKKO
import scipy 
from scipy.optimize import minimize, fsolve

class Solver:
    def __init__(self, model: PowerSystem):
        self.model = model
        self.solved = False
        self.__state_dict = {'solver': None, 'variables': {}}

    def state_dict(self):
        return self.__state_dict
    
    def set_variables(self, vars):
        self.__state_dict['variables'] = vars

    def solve(self, display = True):
        pass

    def extract_results(self):
        if not self.solved:
            raise "Model not solved yet. Call solve()"
        
        variables = self.state_dict()['variables']
        
        # Create dataframes for results
        data = []
        for bus in self.model.buses:
            # Append results as: V, theta, Pload, Qload, Pgen, Qgen
            row = [bus.id, bus.V, bus.theta, bus.Pload, bus.Qload]
            #Pgen = variables['Pgen_fixed'][bus.id][0] + np.sum([variables['Pgen'][gen.bus][0] for gen in self.model.generators if gen.bus == bus.id])
            #Qgen = variables['Qgen_fixed'][bus.id][0] + np.sum([variables['Qgen'][gen.bus][0] for gen in self.model.generators if gen.bus == bus.id])
            Pgen = np.sum([gen.Pgen for gen in self.model.generators if gen.bus == bus.id])
            Qgen = np.sum([gen.Qgen for gen in self.model.generators if gen.bus == bus.id])
            row += [Pgen, Qgen, bus.Vmin, bus.Vmax]

            data.append(row)
        
        # Create DataFrame
        bus_df = pd.DataFrame(data, columns = ['id', 'V', 'theta', 'Pload', 'Qload', 'Pgen', 'Qgen', 'Vmin', 'Vmax'])

        # Create generators DataFrame
        data = []
        for gen in self.model.generators:
            row = [gen.id, variables['lg'][gen.id][0], gen.Pgen, gen.Qgen, gen.cost(gen.Pgen), gen.Pmin, gen.Pmax, gen.Qmin, gen.Qmax]
            data.append(row)
        
        gen_df = pd.DataFrame(data, columns = ['id', 'active', 'Pgen', 'Qgen', 'cost', 'Pmin', 'Pmax', 'Qmin', 'Qmax'])

        #Create lines DataFrame
        data = []
        for line in self.model.lines:
            i = line.from_bus
            k = line.to_bus
            row = [line.id, i, k, variables['l'][line.id][0], variables['Pflow'][i,k][0], variables['Pflow'][k,i][0], variables['Qflow'][i,k][0], variables['Qflow'][k,i][0], line.S(variables['Pflow'][i,k][0], variables['Qflow'][i,k][0]), line.Plosses(variables['Pflow'][i,k][0], variables['Pflow'][k,i][0]), line.mva]
            data.append(row)
        
        lines_df = pd.DataFrame(data, columns = ['id', 'from', 'to', 'active', 'Pout', 'Pin', 'Qout', 'Qin', 'S', 'Ploss', 'mva'])

        return bus_df, gen_df, lines_df
    
    def construct_model_solver(self):
        solver = GEKKO()

        # Voltage variables
        V = solver.Array(solver.Var, dim = (self.model.n_buses,), value = 1.0)

        # Angle variables
        theta = solver.Array(solver.Var, dim = (self.model.n_buses,), lb = -np.pi, ub = np.pi, value = 0.0)

        # Active power generated
        Pgen = solver.Array(solver.Var, dim = (self.model.n_gens,))
        Pgen_fixed = solver.Array(solver.Param, dim = (self.model.n_buses,))
        for bus in self.model.buses:
            Pgen_fixed[bus.id].VALUE = bus.Pgen_fixed

        # Reactive power generated
        Qgen = solver.Array(solver.Var, dim = (self.model.n_gens,))
        Qgen_fixed = solver.Array(solver.Param, dim = (self.model.n_buses,))
        for bus in self.model.buses:
            Qgen_fixed[bus.id].VALUE = bus.Qgen_fixed
        
        # Flow through lines
        Pflow = solver.Array(solver.Var, dim = (self.model.n_buses, self.model.n_buses))
        Qflow = solver.Array(solver.Var, dim = (self.model.n_buses, self.model.n_buses))

        # On/Off for lines
        l = solver.Array(solver.Var, dim = (self.model.n_lines,), lb = 0, ub = 1, value = 1, integer = True)

        # On/Off for generators
        lg = solver.Array(solver.Var, dim = (self.model.n_gens,), lb = 0, ub = 1, value = 1, integer = True)

        state_dict = {
            'solver': solver,
            'variables': {
                'Pgen': Pgen,
                'Qgen': Qgen,
                'Pgen_fixed': Pgen_fixed,
                'Qgen_fixed': Qgen_fixed,
                'V': V,
                'theta': theta,
                'Pflow': Pflow,
                'Qflow': Qflow,
                'l': l,
                'lg': lg
            }
        }

        self.__state_dict = state_dict

        self.construct_optim_constr()
    
    def construct_optim_constr(self):
        assert(self.state_dict()['solver'] is not None, "Construct an optim model first")
        
        m = self.state_dict()['solver']
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
            self.optim_constr_line_P_fromto(line) for line in self.model.lines
        ])
        m.Equations([
            self.optim_constr_line_P_tofrom(line) for line in self.model.lines
        ])

        # # Equations for reactive-power flow through lines
        m.Equations([
            self.optim_constr_line_Q_fromto(line) for line in self.model.lines
        ])
        m.Equations([
            self.optim_constr_line_Q_tofrom(line) for line in self.model.lines
        ])

        # Add constraints for slack bus
        m.Equations([
            self.optim_constr_fixed_bus_angle(bus) for bus in self.model.buses if bus.type == 1
        ])

        m.Equations([
            self.optim_constr_fixed_bus_magnitude(bus) for bus in self.model.buses if bus.type == 1
        ])
    
    def optim_constr_kirchoff_P(self, bus):
        variables = self.state_dict()['variables']
        Pflow = 0.0
        Pgen = variables['Pgen_fixed'][bus.id]
        #Pgen = 0.0

        # Find if there is a generator connected at this bus
        for gen in self.model.generators:
            if gen.bus == bus.id:
                Pgen += variables['Pgen'][gen.id]

        for line in self.model.lines:
            if line.from_bus == bus.id:
                Pflow += variables['Pflow'][line.from_bus, line.to_bus]

        for line in self.model.lines:
            if line.to_bus == bus.id:
                Pflow += variables['Pflow'][line.to_bus, line.from_bus]

        return Pgen == bus.Pload + Pflow
    
    def optim_constr_kirchoff_Q(self, bus):
        variables = self.state_dict()['variables']
        Qflow = 0.0
        Qgen = variables['Qgen_fixed'][bus.id]
        #Qgen = 0.0

        # Find if there is a generator connected at this bus
        for gen in self.model.generators:
            if gen.bus == bus.id:
                Qgen += variables['Qgen'][gen.id]

        for line in self.model.lines:
            if line.from_bus == bus.id:
                Qflow += variables['Qflow'][line.from_bus, line.to_bus]

        for line in self.model.lines:
            if line.to_bus == bus.id:
                Qflow += variables['Qflow'][line.to_bus, line.from_bus]

        return Qgen == bus.Qload + Qflow
    
    def optim_constr_line_P_fromto(self, line):
        m, variables = self.state_dict()['solver'], self.state_dict()['variables']
        l, Pflow, V, theta = variables['l'], variables['Pflow'], variables['V'], variables['theta']
        i = line.from_bus
        k = line.to_bus

        return Pflow[i,k] == l[line.id]*((-self.model.G[i,k] + self.model.g[i,k])*V[i]**2 + V[i]*V[k]*(self.model.G[i,k]*m.cos(theta[i] - theta[k]) + self.model.B[i,k]*m.sin(theta[i] - theta[k])))
    
    def optim_constr_line_P_tofrom(self, line):
        m, variables = self.state_dict()['solver'], self.state_dict()['variables']
        l, Pflow, V, theta = variables['l'], variables['Pflow'], variables['V'], variables['theta']

        i = line.to_bus
        k = line.from_bus

        return Pflow[i,k] == l[line.id]*((-self.model.G[i,k] + self.model.g[i,k])*V[i]**2 + V[i]*V[k]*(self.model.G[i,k]*m.cos(theta[i] - theta[k]) + self.model.B[i,k]*m.sin(theta[i] - theta[k])))
    
    def optim_constr_line_Q_fromto(self, line):
        m, variables = self.state_dict()['solver'], self.state_dict()['variables']
        l, Qflow, V, theta = variables['l'], variables['Qflow'], variables['V'], variables['theta']
        i = line.from_bus
        k = line.to_bus

        return Qflow[i,k] == l[line.id]*((self.model.B[i,k] - self.model.b[i,k])*V[i]**2 + V[i]*V[k]*(-self.model.B[i,k]*m.cos(theta[i] - theta[k]) + self.model.G[i,k]*m.sin(theta[i] - theta[k])))
    
    def optim_constr_line_Q_tofrom(self, line):
        m, variables = self.state_dict()['solver'], self.state_dict()['variables']
        l, Qflow, V, theta = variables['l'], variables['Qflow'], variables['V'], variables['theta']
        i = line.to_bus
        k = line.from_bus

        return Qflow[i,k] == l[line.id]*((self.model.B[i,k] - self.model.b[i,k])*V[i]**2 + V[i]*V[k]*(-self.model.B[i,k]*m.cos(theta[i] - theta[k]) + self.model.G[i,k]*m.sin(theta[i] - theta[k])))

    def optim_constr_fixed_bus_angle(self, bus):
        theta = self.state_dict()['variables']['theta']
        return theta[bus.id] == 0.0
    
    def optim_constr_fixed_bus_magnitude(self, bus):
        V = self.state_dict()['variables']['V']
        return V[bus.id] == bus.V

    def assign_results(self, V, theta, Pgen, Qgen):
        # Assign results
        for bus in self.model.buses:
            bus.V = V[bus.id]
            bus.theta = theta[bus.id]
            bus.Pgen = np.sum([Pgen[gen.id] for gen in self.model.generators if gen.bus == bus.id])
            bus.Qgen = np.sum([Qgen[gen.id] for gen in self.model.generators if gen.bus == bus.id])

        for gen in self.model.generators:
            gen.Pgen = Pgen[gen.id]
            gen.Qgen = Qgen[gen.id]

    def construct_iterative_solver(self, model: PowerSystem, method = "gauss-seidel"):

        if method == "gauss-seidel":
            solver = GaussSeidel(model, IterativeArgs())
        elif method == "scipy" or method == "fsolve":
            solver = ScipyFsolve(model, IterativeArgs())

        self.__state_dict = {
            'solver': solver,
            'variables': {
                'V': np.array([bus.V for bus in model.buses]),
                'theta': np.array([bus.angle for bus in model.buses]),
                'Pgen': np.array([gen.Pgen for gen in model.generators]),
                'Qgen': np.array([gen.Qgen for gen in model.generators])
            }
        }

        V = self.__state_dict['variables']['V']
        theta = self.__state_dict['variables']['theta']

        Pflow = np.zeros((model.n_buses, model.n_buses))
        Qflow = np.zeros_like(Pflow)

        for line in model.lines:
            i = line.from_bus
            k = line.to_bus
            Pflow[i,k] = line.Pflow(np.array([V[i], V[k]]), np.array([theta[i], theta[k]]))
            Qflow[i,k] = line.Qflow(np.array([V[i], V[k]]), np.array([theta[i], theta[k]]))

            i = line.to_bus
            k = line.from_bus
            Pflow[i,k] = line.Pflow(np.array([V[i], V[k]]), np.array([theta[i], theta[k]]))
            Qflow[i,k] = line.Qflow(np.array([V[i], V[k]]), np.array([theta[i], theta[k]]))

        self.__state_dict['variables']['Pflow'] = Pflow
        self.__state_dict['variables']['Qflow'] = Qflow

        return solver
