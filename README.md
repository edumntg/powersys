# powersys
Python package for Power System Analysis, Optimization and Dynamics

# Installation
```pip install powersys```

# Usage
## Building a Power System model using modularity
```python
import powersys as ps
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

print(model)
```

```console
PowerSystem model
Buses: 3        Lines: 3        Generators: 0
Frequency: 60   MVA-Base: 100

Busbar data:
  id    type     V    angle    Pgen    Qgen    Pload    Qload    Vmin    Vmax    Pgen_fixed    Qgen_fixed
----  ------  ----  -------  ------  ------  -------  -------  ------  ------  ------------  ------------
   0       3  1           0     0       0        0        0      0.95    1.05           0             0
   1       2  1.01        0     0.5     0.3      0        0      0.95    1.05           0.5           0.3
   2       0  1           0     0       0        0.5      0.2    0.95    1.05           0             0

Transmission line data:
  id    from_bus    to_bus    R     X  Z              B    a    mva
----  ----------  --------  ---  ----  -----------  ---  ---  -----
   0           0         1  0.1  0.01  (0.1+0.01j)    0    1    inf
   1           0         2  0.1  0.01  (0.1+0.01j)    0    1    inf
   2           1         2  0.1  0.01  (0.1+0.01j)    0    1    inf
```

## Loading Power System IEEE data
