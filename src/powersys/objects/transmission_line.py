import numpy as np
class Line(object):

    PUBLIC_ID = 0

    def __init__(self, id = None, from_bus = None, to_bus = None, R = np.inf, X = np.inf, B = 0.0, a = 1.0, mva = np.inf):
        if id is None:
            id = Line.PUBLIC_ID
            Line.PUBLIC_ID += 1

        self.id = int(id)
        self.from_bus = int(from_bus)
        self.to_bus = int(to_bus)
        self.R = R
        self.X = X
        self.Z = self.R + 1j*self.X
        self.B = B
        self.a = a
        self.mva = mva
        #self.mva = 100

    @property
    def Y(self):
        return 1.0/self.Z
    
    def S(self, P, Q):
        return np.linalg.norm(np.array([P, Q]))
    
    def Plosses(self, Pout, Pin):
        return Pout + Pin
    
    def Pflow(self, V, theta):
        Y = 1/self.Z
        G = np.real(Y)
        B = np.imag(Y)

        return (-G)*V[0]**2 + V[0]*V[1]*(G*np.cos(theta[0] - theta[1]) + B*np.sin(theta[0] - theta[1]))
    
    def Qflow(self, V, theta):
        Y = 1/self.Z
        G = np.real(Y)
        B = np.imag(Y)

        return (B)*V[0]**2 + V[0]*V[1]*(-B*np.cos(theta[0] - theta[1]) + G*np.sin(theta[0] - theta[1]))
    
    def __str__(self):
        return f"{self.id} - [{self.from_bus}, {self.to_bus}, {self.R}, {self.X}, {self.B}, {self.a}, {self.mva}]"
    
    @staticmethod
    def from_dict(dictio):
        return Line(list(dictio.values()))
    
    def as_list(self):
        return list(self.__dict__.values())