import numpy as np
import pandas as pd
from ..objects.object_collection import ObjectCollection
# from ..objects.busbar import Busb
# from ..objects.transmission_line import Line
# from ..objects.synchronous_machine import Generator
from ..objects import *
from ..data.dataloader import DataLoader
from dataclasses import dataclass
from simple_parsing.helpers import Serializable
import json
from typing import Optional
import networkx as nx
import matplotlib.pyplot as plt

@dataclass
class PowerSystemArgs(Serializable):
    buses: Optional[ObjectCollection] = ObjectCollection()
    lines: Optional[ObjectCollection] = ObjectCollection()
    generators: Optional[ObjectCollection] = ObjectCollection()
    f: Optional[int] = 60 # Hz
    mva_base: Optional[float] = 100# in MVA
    v_base: Optional[float] = 100# in kV

class PowerSystem(object):

    SLACK = 3
    PV = 2
    PQ = 0

    def __init__(self,
                 buses = ObjectCollection(),
                 lines = ObjectCollection(),
                 generators = ObjectCollection(),
                 data = None,
                 f = 60,
                 mva_base = 100,
                 v_base = 100):

        self.f = f # Hz        
        self.w = 2*np.pi*self.f
        self.buses = buses
        self.generators = generators
        self.lines = lines

        # If data provided, check that it is an instance of DataLoader
        if data:
            assert isinstance(data, DataLoader), "Invalid type for parameter 'data'. Expected DataLoader."

            self.buses = data.buses
            self.generators = data.generators
            self.lines = data.lines

        self.mva_base = mva_base
        self.v_base = v_base

        self.Ybus = None
        self.Ybus_load = None

        self.b = None
        self.g = None
        self.G = None
        self.B = None

        self.lf_solved = False
        self.opf_solved = False
        self.ta_solved = False

        self.check()
    
    def construct_ybus(self) -> np.array:
        self.Ybus = np.zeros((self.N, self.N), dtype="complex_")
        self.b = np.zeros((self.N, self.N))
        self.g = np.zeros((self.N, self.N))

        # Loop through lines
        for line in self.lines:
            
            self.Ybus[line.from_bus, line.from_bus] += line.Y * (1.0/line.a**2)
            self.Ybus[line.to_bus, line.to_bus] += line.Y * (1.0/line.a**2)

            self.Ybus[line.from_bus, line.to_bus] -= line.Y * 1.0/line.a
            self.Ybus[line.to_bus, line.from_bus] -= line.Y * 1.0/line.a

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
                collection.add(object_class(*row))

            return collection
        except Exception as e:
            raise f"Cannot read CSV {filename}: {e}"

    @staticmethod
    def load_buses(filename):
        if filename.endswith('.csv'): # csvfile
            return PowerSystem.load_from_csv(filename, Busbar)
    
    @staticmethod
    def load_lines(filename):
        if filename.endswith('.csv'):
            return PowerSystem.load_from_csv(filename, Line)
    
    @staticmethod
    def load_gens(filename):
        if filename.endswith('.csv'):
            return PowerSystem.load_from_csv(filename, Generator)
        
    @staticmethod
    def read_json(filename):
        assert filename.endswith('.json'), "Invalid JSON file format."

        with open(filename, 'r') as file:
            data = json.load(file)

            # Read required data
            f = data['f']
            buses = data['buses']
            lines = data['lines']
            generators = data['generators']

            # Create arguments
            args = PowerSystemArgs(
                f = f,
                buses = ObjectCollection().from_dict(buses, Bus),
                lines = ObjectCollection().from_dict(lines, Line),
                generators = ObjectCollection().from_dict(generators, Generator)
            )

            system = PowerSystem(args)
            return system

    def kron_reduction(self):
        # This function calculates the kron reduction of the system
        # The resulting system will be of size (n_gens, n_gens)

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
        if not self.lf_solved:
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

        assert self.lf_solved, "You need to perform a Load Flow analysis first"

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
    def n_pq(self):
        pq_buses = self.buses.filter(lambda x: x.type == PowerSystem.PQ)
        return len(pq_buses)
    
    @property
    def n_pv(self):
        pv_buses = self.buses.filter(lambda x: x.type == PowerSystem.PV)
        return len(pv_buses)

    @property
    def non_PQ_N(self):
        N = 0
        for bus in self.buses:
            if self.get_gen_by_bus(bus.id) or bus.reference: # bus is not PQ
                N += 1

        return N

    @staticmethod
    def from_txt(filename):
        """ NOTE: In order for this function to work, the .txt file must be in IEEE format """

        model = PowerSystem()

        with open(filename, "r") as f:
            # Read all lines
            lines = f.readlines()
        # First line contains date, and others so we pop it
        lines.pop(0)
        # Second line contains title for bus data, so we pop it
        lines.pop(0)
        while True:
            line = lines.pop(0).strip() # current line
            # If line is equal to -999, end loop
            if line == "-999":
                break
                
            # Split
            row = line.split()

            # Get data of interest
            id = int(row[0])
            bus_type = int(row[5])
            V = float(row[6])
            angle = float(row[7])
            Pload = float(row[8])
            Qload = float(row[9])
            Pgen = float(row[10])
            Qgen = float(row[11])

            model.add(Busbar(
                id = id,
                type = bus_type,
                V = V,
                angle = angle,
                Pgen = Pgen,
                Qgen = Qgen,
                Pload = Pload,
                Qload = Qload
            ))

        # Next line contains ile for branch data
        lines.pop(0)
        while True:
            line = lines.pop(0).strip()
            if line == "-999":
                break

            row = line.split()
            from_bus = int(row[0])
            to_bus = int(row[1])
            R = float(row[6])
            X = float(row[7])
            B = float(row[8])

            model.add(Line(
                from_bus = from_bus,
                to_bus = to_bus,
                R = R,
                X = X,
                B = B
            ))

        return model

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
    
    def check(self):
        if len(self.buses):
            return
        
        for bus in self.buses:
            gen = self.get_gen_by_bus(bus.id)
            if gen and bus.type == PowerSystem.PQ:
                print(f"Bus {bus.id} defined as PQ but has generator connected. Switched to PV")
                bus.type = PowerSystem.PV

            if bus.type is None:
                raise Exception(f"Busbar {bus.id} has no valid type. Got: {bus.type}")
        
        return True
    
    def add(self, item):
        """
        Add buses, lines or generators to the model
        """
        if isinstance(item, Busbar):
            self.buses.add(item)
        elif isinstance(item, Line):
            self.lines.add(item)
        elif isinstance(item, Generator):
            self.generators.add(item)
        else:
            raise Exception("Invalid type of Power System object. Got: " + type(item))
        
    def get_graph(self):
        # If no buses or lines, raise exception
        if len(self.buses) == 0 or len(self.lines) == 0:
            raise Exception("No buses or lines specified. Impossible to build graph network")

        # Construct a new graph
        G = nx.Graph()

        # Add buses as nodes
        for bus in self.buses:
            G.add_node(bus.id, P=bus.Pgen, Q=bus.Qgen)

        # Now add edges
        for line in self.lines:
            G.add_edge(line.from_bus, line.to_bus, P = 0.0, Q = 0.0)

        return G
    
    def plot(self, pu = True, zero_index = True, radians = False):
        G = self.get_graph()

        pos = nx.spring_layout(G)

        # Plot
        plt.figure()
        nx.draw(G, with_labels = True)

        angle_multiplier = 1.0
        if radians:
            angle_multiplier = 180.0/np.pi

        mva_multiplier = 1.0
        v_multiplier = 1.0
        if not pu:
            mva_multiplier = self.mva_base
            v_multiplier = self.v_base

        # create labels for each node
        node_labels = {}
        for bus in self.buses:
            node_labels[bus.id] = "V={:.2f}∠{:.2f}\nP={:.4f}\nQ={:.4f}".format(bus.V*v_multiplier, bus.angle*angle_multiplier, bus.Pgen*mva_multiplier, bus.Qgen*mva_multiplier)

        edge_labels = {}
        for line in self.lines:
            from_bus = self.get_bus(line.from_bus)
            to_bus = self.get_bus(line.to_bus)
            P = line.Pflow([from_bus.V, to_bus.V], [from_bus.angle, to_bus.angle])
            Q = line.Qflow([from_bus.V, to_bus.V], [from_bus.angle, to_bus.angle])
            edge_labels[(line.from_bus, line.to_bus)] = "P={:.4f}\nQ={:.4f}".format(P*mva_multiplier, Q*mva_multiplier)

        # node_labels = nx.get_node_attributes(G, 'P')
        nx.draw_networkx_labels(G, pos, labels = node_labels)
        # edge_labels = nx.get_edge_attributes(G, 'P')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        plt.show()

    def __str__(self):
        output = "PowerSystem model\n"
        output += f"Buses: {self.n_buses}\tLines: {self.n_lines}\tGenerators: {self.n_gens}\n"
        output += f"Frequency: {self.f}\tMVA-Base: {self.mva_base}\n"
        output += "\n"
        output += "Busbar data:\n"
        output += str(self.buses)
        output += "\n\n"
        output += "Transmission line data:\n"
        output += str(self.lines)
        output += "\n\n"
        output += "Generators data:\n"
        output += str(self.generators)

        return output

        