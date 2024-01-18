from src.powersys import *
from src.powersys.models import PowerSystem, PowerSystemArgs
from src.powersys.solvers import LFSolver
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = PowerSystemArgs(
    f = 60,
    buses = PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = PowerSystem(args)

# Construct Ybus
Ybus, _, _, _, _ = system.construct_ybus()

# Solve load flow
solver = LFSolver(system)
solver.solve()

# Now, construct the Ybus-load
system.construct_load_ybus()

# Construct Kron
Ykron = system.kron_reduction()

# Construct Yrm
Yrm = system.YRM()

# Construct RM vectors
Vrm, Irm = system.rm()

system.compute_terminal_values()

print([gen.Vt for gen in system.generators])
print([gen.It for gen in system.generators])
print([gen.Ef for gen in system.generators])
print([gen.df for gen in system.generators])


