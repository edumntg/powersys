import sys
sys.path.append('../../src/powersys')

from powersys import power_system, power_system_solver
import pandas as pd

# buses_arr = pd.read_excel('../src/examples/3-bus-system/data/IEEE30.xlsx', 'BUS').to_numpy().astype('float64')
# lines_arr = pd.read_excel('../src/examples/3-bus-system/data/IEEE30.xlsx', 'RAMAS').to_numpy().astype('float64')
# gens_arr = pd.read_excel('../src/examples/3-bus-system/data/IEEE30.xlsx', 'GEN').to_numpy().astype('float64')

# # Preprocess
# buses = dict()
# lines = dict()
# gens = dict()

# for i, row in enumerate(buses_arr):
#    buses[i]= list(row)

# for i, row in enumerate(lines_arr):
#    lines[i] = list(row)

# for i, row in enumerate(gens_arr):
#    gens[i] = list(row)

# data = {
#     'buses': buses,
#     'lines': lines,
#     'generators': gens,
#     'base': 100 # mva
# }
# system = PowerSystem(data)
system = power_system()
system.load_buses('ieee30_buses.csv')
system.load_lines('ieee30_lines.csv')
system.load_gens('ieee30_gens.csv')
# print(str(system))
system.construct_ybus()
print(system.Ybus)
import numpy as np
np.savetxt('ybus.txt', system.Ybus, delimiter=',', fmt='%.4f')
solver = power_system_solver(system)

# Solve
solver.loadflow()
bus, gens, lines = solver.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(bus)
input()
print(gens)
input()
print(lines)


