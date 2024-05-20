import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))
args = ps.models.PowerSystemArgs(
    buses = ps.models.PowerSystem.load_buses(CWD + '/../sample_data/ieee9_buses.csv'),
    lines = ps.models.PowerSystem.load_lines(CWD + '/../sample_data/ieee9_lines.csv'),
    generators = ps.models.PowerSystem.load_gens(CWD + '/../sample_data/ieee9_gens.csv'),
    f = 60 # Hz
)

# system = ps.models.PowerSystem(args)
system = ps.models.PowerSystem(buses = args.buses, lines = args.lines, generators = args.generators)
#system.generators[2].turn_off() # disable generator
print(system.buses)
print(system.lines)
solver = ps.solvers.LF(system)
# Solve
solver.solve()
print(system)