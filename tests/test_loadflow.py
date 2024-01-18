from src.powersys import *
from src.powersys.solvers import LFSolver
from src.powersys.models import PowerSystem, PowerSystemArgs
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = PowerSystemArgs(
    buses = PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv'),
    f = 60 # Hz
)

system = PowerSystem(args)

solver = LFSolver(system)
# Solve
solver.solve()
buses, gens, lines = solver.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(buses)
input()
print(gens)
input()
print(lines)


