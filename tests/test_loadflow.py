from src.powersys import *
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

system = PowerSystem()
system.load_buses(CWD + '/sample_data/ieee30_buses.csv')
system.load_lines(CWD + '/sample_data/ieee30_lines.csv')
system.load_gens(CWD + '/sample_data/ieee30_gens.csv')

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


