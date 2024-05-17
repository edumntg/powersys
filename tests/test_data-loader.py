import powersys as ps
from powersys.data import DataLoader
import os

CWD = os.path.dirname(os.path.abspath(__file__))

ieee9 = ps.datasets.IEEE9()

data = DataLoader(CWD + '/sample_data/ieee9.txt')
data.load()

# Create model
model = ps.models.PowerSystem(data = data)

print(model)