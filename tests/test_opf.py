import sys
sys.path.append('../src/powersys')

from powersys import *
import pandas as pd
import numpy as np

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
system = PowerSystem()
system.load_buses('ieee30_buses.csv')
system.load_lines('ieee30_lines.csv')
system.load_gens('ieee30_gens.csv')
print(str(system))

solver = PowerSystemSolver(system)

# Solve
solver.opf()
bus, gens, lines = solver.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(bus)
input()
print(gens)
input()
print(lines)

print(str(system))

