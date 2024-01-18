from src.powersys import *
from src.powersys.solvers import OPFSolver
from src.powersys.models import PowerSystem, PowerSystemArgs
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = PowerSystemArgs(
    f = 60,
    buses = PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = PowerSystem(args)

solver = OPFSolver(system)

# Solve
solver.solve()
bus, gens, lines = solver.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(bus)
input()
print(gens)
input()
print(lines)

print(str(system))

