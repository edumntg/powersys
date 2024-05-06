import numpy as np

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

        # Generator values like impedances, reactances, etc
        self.Ra = data[9]
        self.Xd = data[10]
        self.Xq = data[11]
        self.Xdp = data[12]
        self.Xqp = data[13]
        self.Xdpp = data[14]
        self.Xqpp = data[15]
        self.Td0p = data[16]
        self.Tq0p = data[17]
        self.Td0pp = data[18]
        self.Tq0pp = data[19]
        self.H = data[20]
        self.USE_CAGE = data[21] == 1
        self.USE_AGC = data[22] == 1
        self.Kv = data[23]
        self.Ki = data[24]
        self.Kt = data[25]
        self.R = data[26]
        self.Tv = data[27]
        self.Tt = data[28]
        self.USE_AVR = data[29] == 1
        self.Kmed = data[30]
        self.Kexc = data[31]
        self.Ka = data[32]
        self.Tmed = data[33]
        self.Texc = data[34]
        self.Ta = data[35]
        self.Kd = data[36]
        self.Kp = data[37]
        self.Kvi = data[38]
        self.Vexc_min = data[39]
        self.Vexc_max = data[40]
        self.USE_PSS = data[41] == 1
        self.Kest = data[42]
        self.Tw = data[43]
        self.T1 = data[44]
        self.T2 = data[45]
        self.T3 = data[46]
        self.T4 = data[47]
        self.Vpc2_min = data[48]
        self.Vpc2_max = data[49]

        self.Vt = 0
        self.It = 0
        self.Ef = 0
        self.df = 0
        self.Pm = 0

        self.Vq = 0
        self.Vd = 0
        self.Iq = 0
        self.Id = 0

        self.Eq = 0
        self.Eqp = 0
        self.Eqpp = 0

        self.Ed = 0
        self.Edp = 0
        self.Edpp = 0
        
        self.Eexc = 0

        self.active = 1

    def turn_off(self):
        self.active = 0

    def turn_on(self):
        self.active = 1

    def compute_terminal_values(self, Vt, It):
        self.Vt = Vt
        self.It = It

        self.Ef = self.Vt + (self.Ra + 1j*self.Xq)*self.It
        self.df = np.angle(self.Ef)

    def compute_gap_power(self):
        self.Pm = np.real(self.Vt*np.conj(self.It)) + self.Ra*np.abs(self.It)**2.0

    def compute_excitation_values(self):
        self.Eq = self.Vq + self.Ra*self.Iq - self.Xd*self.Id
        self.Eqp = self.Vq + self.Ra*self.Iq - self.Xdp*self.Id
        self.Eqpp = self.Vq + self.Ra*self.Iq - self.Xdpp*self.Id

        self.Ed = self.Vd + self.Ra*self.Id + self.Xq*self.Iq
        self.Edp = self.Vd + self.Ra*self.Id + self.Xqp*self.Iq
        self.Edpp = self.Vd + self.Ra*self.Id + self.Xqpp*self.Iq

        self.Eexc = np.abs(self.Eq + 1j*self.Ed)

        self.Vvi = self.Eexc / (self.Kexc*self.Ka)
        self.Va = self.Eexc / self.Kexc

    def cost(self, P):
        return self.a*P**2 + self.b*P + self.c
    
    def __str__(self):
        return f"{self.id} - [{self.bus}, {self.Pgen:.2f}, {self.Qgen:.2f}, {self.c}, {self.b}, {self.a}, {self.Pmin:.2f}, {self.Pmax:.2f}, {self.Qmin:.2f}, {self.Qmax:.2f}]"
    
    @staticmethod
    def from_dict(dictio):
        return Generator(list(dictio.values()))