import powersys as ps
import os

CWD = os.path.dirname(os.path.abspath(__file__))

model = ps.models.PowerSystem.from_txt(CWD + '/sample_data/ieee9.txt')
model.plot()