import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = ps.models.PowerSystemArgs(
    buses = ps.models.PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = ps.models.PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = ps.models.PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv'),
    f = 60 # Hz
)

system = ps.models.PowerSystem(args)

#system.generators[2].turn_off() # disable generator

solver = ps.solvers.LF(system)
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