import numpy as np
import pandas as pd
from ..objects.ObjectCollection import ObjectCollection
from ..objects.Bus import Bus
from ..objects.Line import Line
from ..objects.Generator import Generator
from dataclasses import dataclass
from simple_parsing.helpers import Serializable

@dataclass
class PowerSystemArgs(Serializable):
    buses: ObjectCollection
    lines: ObjectCollection
    generators: ObjectCollection
    f: int = 60 # Hz


class PowerSystem(object):

    def __init__(self, args: PowerSystemArgs):
        self.args = args

        self.f = self.args.f # Hz
        self.w = 2*np.pi*self.f
        self.buses = args.buses
        self.generators = args.generators
        self.lines = args.lines

        self.Ybus = None
        self.Ybus_load = None

        self.b = None
        self.g = None
        self.G = None
        self.B = None

        self.lf_solved = False
        self.opf_solved = False
        self.ta_solved = False
    
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

    @staticmethod
    def load_from_csv(filename, object_class, delimiter = ','):
        # Open csv file
        try:
            data = pd.read_csv(filename, delimiter = delimiter).to_numpy().astype(float)
            # Take each row and create a new object_class element
            collection = ObjectCollection()
            for row in data:
                collection.add(object_class(row))

            return collection
        except Exception as e:
            raise f"Cannot read CSV {filename}: {e}"

    @staticmethod
    def load_buses(filename):
        if filename.endswith('.csv'): # csvfile
            return PowerSystem.load_from_csv(filename, Bus)
    
    @staticmethod
    def load_lines(filename):
        if filename.endswith('.csv'):
            return PowerSystem.load_from_csv(filename, Line)
    
    @staticmethod
    def load_gens(filename):
        if filename.endswith('.csv'):
            return PowerSystem.load_from_csv(filename, Generator)

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

    def compute_terminal_values(self):
        # Calculate the voltages and currents at bus terminals

        if self.Vrm is None or self.Irm is None:
            raise "No RM values calculated"
        
        for i, gen in enumerate(self.generators):
            Vt = self.Vrm[2*i] +1j*self.Vrm[2*i+1]
            It = self.Irm[2*i] + 1j*self.Irm[2*i+1]

            gen.compute_terminal_values(Vt[0], It[0])

    def compute_gap_power(self):
        # Calculates the mechanical power supplied to the generator (electrical + losses)

        for i, gen in enumerate(self.generators):
            gen.compute_gap_power()

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
                    T[2*i, 2*i] = np.cos(geni.df)
                    T[2*i+1, 2*i+1] = np.cos(geni.df)

                    T[2*i, 2*i+1] = np.sin(geni.df)
                    T[2*i+1, 2*i] = -np.sin(geni.df)

        self.T = T

        return self.T
    
    def park_transform(self, M):
        if M is None:
            M = self.m_reduction()

        T = self.park_matrix()
        
        A = np.linalg.inv(np.linalg.inv(T)@M)@T

        return A
    
    def compute_qd_values(self):
        # This function calculates voltages and currents in QD domain
        if self.Vrm is None:
            raise "You need to perform an RM transform first"

        T = self.park_matrix()

        # Calculate Vqd and Iqd
        Vqd = T@self.Vrm
        Iqd = T@self.Irm

        for i, gen in enumerate(self.generators):
            gen.Vq = Vqd[2*i][0]
            gen.Vd = Vqd[2*i+1][0]
            gen.Iq = Iqd[2*i][0]
            gen.Id = Iqd[2*i+1][0]

    def compute_excitation_values(self):
        # Calculates the internal generator excitation voltages in QD domain

        for i, gen in enumerate(self.generators):
            gen.compute_excitation_values()

    def transient_initial_condition(self):
        # Construct the vectors with initial values for transient analysis
        self.w0 = np.zeros((self.n_gens, 1))
        self.d0 = np.array([gen.df for gen in self.generators])
        self.Pe_gap0 = np.array([gen.Pm for gen in self.generators])
        self.Vmed0 = np.array([np.sqrt(gen.Vq**2.0 + gen.Vd**2.0) for gen in self.generators])
        self.Vexc0 = np.array([gen.Eexc for gen in self.generators])
        self.Vt0 = np.array([bus.V for bus in self.buses])
        self.Vpc10 = np.zeros((self.n_gens, 1))
        self.Vpc20 = np.zeros((self.n_gens, 1))

        self.Vref = self.Vt0.copy()

    def prepare_transient_analysis(self):

        assert self.load_flow_solved, "You need to perform a Load Flow analysis first"

        """ Compute all system's required parameters for transient analysis """
        # Load-ybus
        self.construct_load_ybus()

        # Kron values
        self.kron_reduction()

        # RM Ybus
        self.YRM()

        # RM values
        self.rm()

        # VI terminal values
        self.compute_terminal_values()

        # Gap power
        self.compute_gap_power()

        # M reduction
        self.m_reduction()

        # Park matrix
        self.park_matrix()

        # Qd values
        self.compute_qd_values()

        # Excitation values
        self.compute_excitation_values()

        # Transient initial condition
        self.transient_initial_condition()


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
    
