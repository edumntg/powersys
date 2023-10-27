class Bus(object):

    def __init__(self, data):
        self.id = int(data[0])
        self.reference = int(data[1]) == 1
        self.V = data[2]
        self.theta = data[3]
        self.Pgen = data[4]
        self.Qgen = data[5]
        self.Pload = data[6]
        self.Qload = data[7]
        self.Vmin = data[8]
        self.Vmax = data[9]

        # Fixed generation
        self.Pgen_fixed = data[4]
        self.Qgen_fixed = data[5]

    def __str__(self):
        return f"{self.id} - [{self.V:.2f}V, {self.theta:.2f}rads, P={self.Pgen:.2f}, Q={self.Qgen:.2f}, Vmin={self.Vmin:.2f}, Vmax={self.Vmax:.2f}]"

    