from src.powersys import *
from src.powersys.models import PowerSystem
from src.powersys.solvers import LFSolver
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

system = PowerSystem()
system.load_buses(CWD + '/sample_data/ieee9_buses.csv')
system.load_lines(CWD + '/sample_data/ieee9_lines.csv')
system.load_gens(CWD + '/sample_data/ieee9_gens.csv')

# Construct Ybus
Ybus, _, _, _, _ = system.construct_ybus()

# Solve load flow
solver = LFSolver(system)
solver.solve()

# Now, construct the Ybus-load
system.construct_load_ybus()

# Construct Kron
Ykron = system.kron_reduction()
print(Ykron)



