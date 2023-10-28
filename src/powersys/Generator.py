class Generator(object):
    def __init__(self, data):
        self.id = int(data[0])
        self.bus = int(data[1])
        self.c = data[2]
        self.b = data[3]
        self.a = data[4]
        self.Pmin = data[5]
        self.Pmax = data[6]
        self.Qmin = data[7]
        self.Qmax = data[8]
        self.Pgen = 0
        self.Qgen = 0

    def cost(self, P):
        return self.a*P**2 + self.b*P + self.c
    
    def __str__(self):
        return f"{self.id} - [{self.bus}, {self.Pgen:.2f}, {self.Qgen:.2f}, {self.c}, {self.b}, {self.a}, {self.Pmin:.2f}, {self.Pmax:.2f}, {self.Qmin:.2f}, {self.Qmax:.2f}]"