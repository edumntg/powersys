import powersys as ps

ieee9 = ps.datasets.IEEE9()

model = ps.models.PowerSystem(data = ieee9)
print(model)