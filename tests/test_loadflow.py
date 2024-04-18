import powersys as ps
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = ps.model.PowerSystemArgs(
    buses = ps.model.PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = ps.model.PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = ps.model.PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv'),
    f = 60 # Hz
)

system = ps.model.PowerSystem(args)

solvers = ps.solver.LF(system)
# Solve
solvers.solve()
buses, gens, lines = solvers.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(buses)
input()
print(gens)
input()
print(lines)


