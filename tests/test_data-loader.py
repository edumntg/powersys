import powersys as ps
from powersys.data import DataLoader
import os

CWD = os.path.dirname(os.path.abspath(__file__))

data = DataLoader(CWD + '/sample_data/ieee9.txt')
data.load()

# Create model
model = ps.models.PowerSystem(data = data)

# Create a LF solver
solver = ps.solvers.LF(model)
solver.solve(method="newton-raphson")

model.plot()