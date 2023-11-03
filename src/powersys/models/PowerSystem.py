import numpy as np
import pandas as pd
from ..objects.ObjectCollection import ObjectCollection
from ..objects.Bus import Bus
from ..objects.Line import Line
from ..objects.Generator import Generator

class PowerSystem(object):

    def __init__(self, data: dict = {}):
        self.f = 60 # Hz
        self.we = 2*np.pi*self.f
        self.buses_raw = ObjectCollection()
        self.lines_raw = ObjectCollection()
        self.generators_raw = ObjectCollection()
        self.buses = ObjectCollection()
        self.lines = ObjectCollection()
        self.generators = ObjectCollection()

        if data and type(data) is dict and len(data) > 0:
            self.buses_raw = data['buses']
            self.lines_raw = data['lines']
            self.generators_raw = data['generators']

            # Initialize buses, lines and generators
            for id, data in self.buses_raw.items():
                bus = Bus(data)
                self.buses.append(bus)

            for id, data in self.lines_raw.items():
                line = Line([id] + data)
                self.lines.append(line)

            for id, data in self.generators_raw.items():
                gen = Generator([id] + data)
                self.generators.append(gen)

        self.Ybus = None

        self.b = None
        self.g = None
        self.G = None
        self.B = None

        self.load_flow_solved = False
        self.Ybus_load = None

    def construct_ybus(self) -> np.array:
        self.Ybus = np.zeros((self.N, self.N), dtype="complex_")
        self.b = np.zeros((self.N, self.N))
        self.g = np.zeros((self.N, self.N))

        # Loop through lines
        for line in self.lines:

            self.Ybus[line.from_bus, line.from_bus] += line.Y() * (1.0/line.a**2)
            self.Ybus[line.to_bus, line.to_bus] += line.Y() * (1.0/line.a**2)

            self.Ybus[line.from_bus, line.to_bus] -= line.Y() * 1.0/line.a
            self.Ybus[line.to_bus, line.from_bus] -= line.Y() * 1.0/line.a

            # Add shunts
            self.Ybus[line.from_bus, line.from_bus] += 1j*line.B/2.0
            self.Ybus[line.to_bus, line.to_bus] += 1j*line.B/2.0

            self.b[line.from_bus, line.to_bus] = line.B/2.0
            self.b[line.to_bus, line.from_bus] = line.B/2.0
            
        self.G = self.Ybus.real
        self.B = self.Ybus.imag

        return self.Ybus, self.G, self.B, self.g, self.b
    
    def get_gen(self, id):
        for gen in self.generators:
            if gen.id == id:
                return gen
            
        return None
    
    def get_gen_by_bus(self, bus_id):
        for gen in self.generators:
            if gen.bus == bus_id:
                return gen
        
        return None

    def get_bus(self, id):
        for bus in self.buses:
            if bus.id == id:
                return bus
        
        return None

    def get_line(self, id):
        for line in self.lines:
            if line.id == id:
                return line
            
        return None

    def load_from_csv(self, filename, object_class, delimiter = ','):
        # Open csv file
        try:
            data = pd.read_csv(filename, delimiter = delimiter).to_numpy().astype(float)
            # Take each row and create a new object_class element
            collection = ObjectCollection()
            for row in data:
                collection.add(object_class(row))

            if object_class is Bus:
                self.buses = collection
            elif object_class is Line:
                self.lines = collection
            elif object_class is Generator:
                self.generators = collection
        except Exception as e:
            raise f"Cannot read CSV {filename}: {e}"

    def load_buses(self, filename):
        if filename.endswith('.csv'): # csvfile
            return self.load_from_csv(filename, Bus)
        
    def load_lines(self, filename):
        if filename.endswith('.csv'):
            return self.load_from_csv(filename, Line)
    
    def load_gens(self, filename):
        if filename.endswith('.csv'):
            return self.load_from_csv(filename, Generator)

    def kron_reduction(self):
        # This function calculates the cron reduction of the system

        # Ya is a matrix of size (n_gens, n_gens)
        Ya = self.Ybus_load[:self.n_gens, :self.n_gens]

        # Yb is a matrix of size (n_gens, n_buses-n_gens)
        Yb = self.Ybus_load[:self.n_gens, self.n_gens:]

        # Yc is the transpose of Yb
        Yc = Yb.T

        # Yd is a matrix that goes from n_gens: and n_gens:
        Yd = self.Ybus_load[self.n_gens:, self.n_gens:]

        # Now, calculate Kron eduction
        Ykron = Ya - Yb@(np.linalg.inv(Yd)@Yc)

        self.Ykron = Ykron

        return Ykron

    def construct_load_ybus(self):
        # This function creates a Ybus and adds bus loads as impedances
        if not self.load_flow_solved:
            raise "You need to perform a Load Flow study first"

        # Now, calculate load impedances at each bus
        self.Ybus_load = self.Ybus.copy()
        for bus in self.buses:
            S = bus.Pload + 1j*bus.Qload;
            if np.linalg.norm(S) > 0.0:
                Z = bus.V**2 / np.conj(S)

                self.Ybus_load[bus.id,bus.id] += 1.0 / Z;
    
        return self.Ybus_load

    def YRM(self):
        # Create the RM equivalent matrix
        Gk = np.real(self.Ykron)
        Bk = np.imag(self.Ykron)

        nk = self.Ykron.shape[0]

        Yrm = np.zeros((2*nk, 2*nk))

        for i in range(nk):
            for j in range(nk):
                if i == j:
                    Yrm[2*i,2*i] = Gk[i,i]
                    Yrm[2*i+1,2*i+1] = Gk[i,i]

                    Yrm[2*i,2*i+1] = -Bk[i,i]
                    Yrm[2*i+1, 2*i] = Bk[i,i]
                else:
                    Yrm[2*i,2*j] = Gk[i,j]
                    Yrm[2*i+1,2*j+1] = Gk[i,j]

                    Yrm[2*i,2*j+1] = -Bk[i,j]
                    Yrm[2*i+1,2*j] = Bk[i,j]

        self.Yrm = Yrm
        return self.Yrm
    
    def rm(self):
        # Calculates the Vrm and Irm vectors
        if self.Yrm is None:
            raise "No YRM matrix created"
        
        Vrm = np.zeros((2*self.non_PQ_N, 1))
        Irm = np.zeros((2*self.non_PQ_N, 1))

        for bus in self.buses:
            if self.get_gen_by_bus(bus.id) or bus.reference: # not PQ
                # Calculate complex expression for voltage at this bus
                V_rectangular = bus.V*(np.cos(bus.theta) + 1j*np.sin(bus.theta))
                Vrm[2*bus.id] = np.real(V_rectangular)
                Vrm[2*bus.id+1] = np.imag(V_rectangular)

        Irm = self.Yrm@Vrm

        self.Vrm = Vrm
        self.Irm = Irm
        
        return Vrm, Irm

    def vi_terminal_values(self):
        # Calculate the voltages and currents at bus terminals

        if self.Vrm is None or self.Irm is None:
            raise "No RM values calculated"
        
        Vt = np.zeros((self.n_gens, 1), dtype='complex_')
        It = np.zeros((self.n_gens, 1), dtype='complex_')

        # Vector to store internal generator voltages
        Ef = np.zeros((self.n_gens, 1), dtype='complex_')
        df = np.zeros((self.n_gens, 1))
        for i, gen in enumerate(self.generators):
            Vt[i] = self.Vrm[2*i] +1j*self.Vrm[2*i+1]
            It[i] = self.Irm[2*i] + 1j*self.Irm[2*i+1]

            Ef[i] = Vt[i] + (gen.Ra + 1j*gen.Xq)*It[i]
            df[i] = np.angle(Ef[i])

            gen.Vt = Vt[i]
            gen.It = It[i]

            gen.Ef = Ef[i]
            gen.df = df[i]

        self.Vt = Vt
        self.It = It
        self.Ef = Ef
        self.df = df

        return self.Vt, self.It, self.Ef, self.df

    def gap_power(self):
        # Calculates the mechanical power supplied to the generator (electrical + losses)

        Pm = np.zeros((self.n_gens, 1))
        for i, gen in enumerate(self.generators):
            Pm[i] = np.real(self.Vt[i]*np.conj(self.It[i])) + gen.Ra*np.abs(self.It[i])**2.0
            gen.Pm = Pm[i]
        
        self.Pm = Pm

        return self.Pm

    def m_reduction(self):
        # Reduces the system by removing PQ buses

        M = np.zeros((2*self.n_gens, 2*self.n_gens))

        for i, gen in enumerate(self.generators):
            for j, gen in enumerate(self.generators):
                if i == j:
                    M[2*i, 2*i] = gen.Ra
                    M[2*i+1, 2*i+1] = gen.Ra

                    Xd = gen.Xdp
                    Xq = gen.Xqp

                    M[2*i,2*i+1] = -Xd
                    M[2*i+1,2*i] = Xq

        self.M_reduction = M

        return self.M_reduction

    def park_matrix(self):
        if self.M_reduction is None:
            raise "You need to perform an M reduction first"
        
        T = np.zeros((2*self.n_gens, 2*self.n_gens))
        for i, geni in enumerate(self.generators):
            for j, genj in enumerate(self.generators):
                if i == j:
                    T[2*i, 2*i] = np.cos(self.df[geni.bus])
                    T[2*i+1, 2*i+1] = np.cos(self.df[geni.bus])

                    T[2*i, 2*i+1] = np.sin(self.df[geni.bus])
                    T[2*i+1, 2*i] = -np.sin(self.df[geni.bus])

        self.T = T

        return self.T
    
    def park_transform(self, M):
        if M is None:
            M = self.m_reduction()

        T = self.park_matrix()
        
        A = np.linalg.inv(np.linalg.inv(T)@M)@T

        return A
    
    def qd_values(self):
        # This function calculates voltages and currents in QD domain
        if self.Vrm is None:
            raise "You need to perform an RM transform first"

        T = self.park_matrix()

        # Calculate Vqd and Iqd
        Vqd = T@self.Vrm
        Iqd = T@self.Irm

        Vq = np.zeros((self.n_gens, 1))
        Vd = np.zeros((self.n_gens, 1))
        Iq = np.zeros((self.n_gens, 1))
        Id = np.zeros((self.n_gens, 1))

        for i, gen in enumerate(self.generators):
            Vq[i] = Vqd[2*i]
            Vd[i] = Vqd[2*i+1]
            Iq[i] = Iqd[2*i]
            Id[i] = Iqd[2*i+1]

            gen.Vq = Vq[i]
            gen.Vd = Vd[i]
            gen.Iq = Iq[i]
            gen.Id = Id[i]

        self.Vq = Vq
        self.Vd = Vd
        self.Iq = Iq
        self.Id = Id

        return Vq, Vd, Iq, Id

    def excitation_values(self):
        # Calculates the internal generator excitation voltages in QD domain
        Eq = np.zeros((self.n_gens, 1))
        Eqp = np.zeros((self.n_gens, 1))
        Eqpp = np.zeros((self.n_gens, 1))

        Ed = np.zeros((self.n_gens, 1))
        Edp = np.zeros((self.n_gens, 1))
        Edpp = np.zeros((self.n_gens, 1))

        Eexc = np.zeros((self.n_gens, 1), dtype = "complex_")
        for i, gen in enumerate(self.generators):
            Eq[i] = self.Vq[i] + gen.Ra*self.Iq[i] - gen.Xd*self.Id[i]
            Eqp[i] = self.Vq[i] + gen.Ra*self.Iq[i] - gen.Xdp*self.Id[i]
            Eqpp[i] = self.Vq[i] + gen.Ra*self.Iq[i] - gen.Xdpp*self.Id[i]

            Ed[i] = self.Vd[i] + gen.Ra*self.Id[i] + gen.Xq*self.Iq[i]
            Edp[i] = self.Vd[i] + gen.Ra*self.Id[i] + gen.Xqp*self.Iq[i]
            Edpp[i] = self.Vd[i] + gen.Ra*self.Id[i] + gen.Xqpp*self.Iq[i]

            Eexc[i] = np.abs(Eq[i] + 1j*Ed[i])

            gen.Eq = Eq[i]
            gen.Eqp = Eqp[i]
            gen.Eqpp = Eqpp[i]

            gen.Ed = Ed[i]
            gen.Edp = Edp[i]
            gen.Edpp = Edpp[i]

            gen.Eexc = Eexc[i]
            gen.Vvi = gen.Eexc / (gen.Kexc*gen.Ka)
            gen.Va = gen.Eexc / gen.Kexc

        self.Eq = Eq
        self.Eqp = Eqp
        self.Eqpp = Eqpp

        self.Ed = Ed
        self.Edp = Edp
        self.Edpp = Edpp

        self.Eexc = Eexc

        return Eq, Ed, Eexc

    def transient_initial_condition(self):
        # Construct the vectors with initial values for transient analysis
        self.w0 = np.zeros((self.n_gens, 1))
        self.d0 = self.df.copy()
        self.Pe_gap0 = self.Pm.copy()
        self.Vmed0 = np.sqrt(self.Vq**2.0 + self.Vd**2.0).copy()
        self.Vexc0 = self.Eexc.copy()
        self.Vt0 = self.V.copy()
        self.Vpc10 = np.zeros((self.n_gens, 1))
        self.Vpc20 = np.zeros((self.n_gens, 1))

        self.Vref = self.V.copy()

    def transient_analysis(self):
        pass


    @property
    def N(self):
        return len(self.buses)

    @property
    def n_gens(self):
        return len(self.generators)
    
    @property
    def n_buses(self):
        return len(self.buses)
    
    @property
    def n_lines(self):
        return len(self.lines)

    @property
    def non_PQ_N(self):
        N = 0
        for bus in self.buses:
            if self.get_gen_by_bus(bus.id) or bus.reference: # bus is not PQ
                N += 1

        return N

    @staticmethod
    def from_excel(filename):
        pass

    @staticmethod
    def from_pandas(df):
        pass

    @staticmethod
    def from_dir(df):
        pass

    @staticmethod
    def from_json(df):
        pass
        
    def __str__(self):
        return f"Buses:\n{str(self.buses)}\n\nLines:\n{str(self.lines)}\n\nGenerators:\n{str(self.generators)}"
    
