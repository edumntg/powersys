from powersys.powersys import *
import pandas as pd
import numpy as np

system = PowerSystem()
system.load_buses('ieee30_buses.csv')
system.load_lines('ieee30_lines.csv')
system.load_gens('ieee30_gens.csv')

solver = PowerSystemSolver(system)
# Solve
solver.loadflow()
buses, gens, lines = solver.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(buses)
input()
print(gens)
input()
print(lines)


