import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = ps.model.PowerSystemArgs(
    f = 60,
    buses = ps.model.PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = ps.model.PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = ps.model.PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = ps.model.PowerSystem(args)

solvers = ps.solver.OPF(system)

# Solve
solvers.solve()
bus, gens, lines = solvers.extract_results()
print("Optimization problem solved. Press enter to continue...")
input()
print(bus)
input()
print(gens)
input()
print(lines)

print(str(system))

