import numpy as np

class Generator(object):
    PUBLIC_ID = 0
    def __init__(self, 
             id, 
             bus, 
             c, 
             b, 
             a, 
             Pmin, 
             Pmax, 
             Qmin, 
             Qmax, 
             Ra, 
             Xd, 
             Xq, 
             Xdp, 
             Xqp, 
             Xdpp, 
             Xqpp, 
             Td0p, 
             Tq0p, 
             Td0pp, 
             Tq0pp, 
             H, 
             USE_CAGE, 
             USE_AGC, 
             Kv, 
             Ki, 
             Kt, 
             R, 
             Tv, 
             Tt, 
             USE_AVR, 
             Kmed, 
             Kexc, 
             Ka, 
             Tmed, 
             Texc, 
             Ta, 
             Kd, 
             Kp, 
             Kvi, 
             Vexc_min, 
             Vexc_max, 
             USE_PSS, 
             Kest, 
             Tw, 
             T1, 
             T2, 
             T3, 
             T4, 
             Vpc2_min, 
             Vpc2_max):
        
        if id is None:
            id = Generator.PUBLIC_ID
            Generator.PUBLIC_ID += 1
                
        self.id = int(id)
        self.bus = int(bus)
        self.c = c
        self.b = b
        self.a = a
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Qmin = Qmin
        self.Qmax = Qmax
        self.Pgen = 0
        self.Qgen = 0

        # Generator values like impedances, reactances, etc
        self.Ra = Ra
        self.Xd = Xd
        self.Xq = Xq
        self.Xdp = Xdp
        self.Xqp = Xqp
        self.Xdpp = Xdpp
        self.Xqpp = Xqpp
        self.Td0p = Td0p
        self.Tq0p = Tq0p
        self.Td0pp = Td0pp
        self.Tq0pp = Tq0pp
        self.H = H
        self.USE_CAGE = USE_CAGE == 1
        self.USE_AGC = USE_AGC == 1
        self.Kv = Kv
        self.Ki = Ki
        self.Kt = Kt
        self.R = R
        self.Tv = Tv
        self.Tt = Tt
        self.USE_AVR = USE_AVR == 1
        self.Kmed = Kmed
        self.Kexc = Kexc
        self.Ka = Ka
        self.Tmed = Tmed
        self.Texc = Texc
        self.Ta = Ta
        self.Kd = Kd
        self.Kp = Kp
        self.Kvi = Kvi
        self.Vexc_min = Vexc_min
        self.Vexc_max = Vexc_max
        self.USE_PSS = USE_PSS == 1
        self.Kest = Kest
        self.Tw = Tw
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.Vpc2_min = Vpc2_min
        self.Vpc2_max = Vpc2_max

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
    
    def as_list(self):
        return list(self.__dict__.values())