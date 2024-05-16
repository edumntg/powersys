import math
class Busbar(object):

    PUBLIC_ID = 0

    def __init__(
            self,
            id = None,
            type = None,
            V = 1.0,
            angle = 0.0,
            Pgen = 0.0,
            Qgen = 0.0,
            Pload = 0.0,
            Qload = 0.0,
            Vmin = 0.95,
            Vmax = 1.05
    ):
        
        if id is None:
            id = Busbar.PUBLIC_ID
            Busbar.PUBLIC_ID += 1
        
        self.id = int(id)
        self.type = type
        self.V = V
        self.angle = angle
        self.Pgen = Pgen
        self.Qgen = Qgen
        self.Pload = Pload
        self.Qload = Qload
        self.Vmin = Vmin
        self.Vmax = Vmax

        self.Pgen_fixed = Pgen
        self.Qgen_fixed = Qgen
    
    @property
    def P(self):
        return self.Pgen

    @property
    def Q(self):
        return self.Qgen

    @property
    def voltage(self):
        return self.V

    @property
    def theta(self):
        return self.angle
    
    @property
    def S(self):
        return math.sqrt(self.P**2.0 + self.Q**2.0)
    
    @property
    def load(self):
        return self.Pload + 1j*self.Qload
    
    @staticmethod
    def from_dict(dictio):
        return Busbar(list(dictio.values()))
    
    def as_list(self):
        return list(self.__dict__.values())

    def __str__(self):
        return f"{self.id} - [{self.V:.2f}V, {self.theta:.2f}rads, P={self.Pgen:.2f}, Q={self.Qgen:.2f}, Vmin={self.Vmin:.2f}, Vmax={self.Vmax:.2f}]"

    