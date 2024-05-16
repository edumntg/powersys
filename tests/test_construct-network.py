import powersys as ps
import os
from powersys.objects import *

# Construct new model
model = ps.models.PowerSystem()

# Add 3 buses
model.add(Busbar(type = model.SLACK, V = 1.00, angle = 0.0)) # Slack bus
model.add(Busbar(type = model.PV, V = 1.01, angle = 0.0, Pgen = 0.5, Qgen = 0.3)) # PV
model.add(Busbar(type = model.PQ, V = 1.00, angle = 0.0, Pload = 0.5, Qload = 0.2)) # PQ

# Add 3 lines to form a ring
model.add(Line(from_bus = 0, to_bus = 1, R = 0.1, X = 0.01))
model.add(Line(from_bus = 0, to_bus = 2, R = 0.1, X = 0.01))
model.add(Line(from_bus = 1, to_bus = 2, R = 0.1, X = 0.01))

# Solve model
solver = ps.solver.LF(model)
solver.solve(method = 'newton-raphson', tol = 1e-6)

model.plot()