from src.powersys.solvers import TransientAnalysisSolver, LF
from src.powersys.models import PowerSystem, PowerSystemArgs
import os
import numpy as np
import matplotlib.pyplot as plt

CWD = os.path.dirname(os.path.abspath(__file__))

args = PowerSystemArgs(
    f = 60,
    buses = PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = PowerSystem(args)

# Construct Ybus
system.construct_ybus()

# Solve load flow
solver = LF(system)
solver.solve()

# Create a new transient stability solver
transient_solver = TransientAnalysisSolver(system)
transient_solver.tspan = np.linspace(0, 50/system.f, 1000) # cycles

# Solve
transient_solver.solve()

plt.figure()
plt.plot(transient_solver.omega[2,:])
plt.show()