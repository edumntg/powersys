import powersys as ps
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = ps.model.PowerSystemArgs(
    f = 60,
    buses = ps.model.PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = ps.model.PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = ps.model.PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = ps.model.PowerSystem(args)

# Construct Ybus
Ybus, _, _, _, _ = system.construct_ybus()

# Solve load flow
solvers = ps.solver.LF(system)
solvers.solve()

# Now, construct the Ybus-load
system.construct_load_ybus()

# Construct Kron
Ykron = system.kron_reduction()

# Construct Yrm
Yrm = system.YRM()

# Construct RM vectors
Vrm, Irm = system.rm()

system.compute_terminal_values()

system.compute_gap_power()

M = system.m_reduction()

T = system.park_matrix()

print(T)



