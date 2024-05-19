import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))

args = ps.models.PowerSystemArgs(
    f = 60,
    buses = ps.models.PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
    lines = ps.models.PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
    generators = ps.models.PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
)

system = ps.models.PowerSystem(buses= args.buses, lines = args.lines, generators = args.generators)

solvers = ps.solvers.OPF(system)

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

system.plot(pu=False, radians=False)