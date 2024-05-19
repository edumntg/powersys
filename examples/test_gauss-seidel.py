import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))

data = ps.datasets.IEEE9()

model = ps.models.PowerSystem(data = data)

print(model)
solver = ps.solvers.LF(model)

# Solve
solver.solve(method = "gauss-seidel", tol=1E-9)

print(model)