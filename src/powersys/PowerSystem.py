import numpy as np
import pandas as pd
from .ObjectCollection import ObjectCollection
from .Bus import Bus
from .Line import Line
from .Generator import Generator

class PowerSystem(object):

    def __init__(self, data: dict = {}):
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

        self.n_buses = len(self.buses)
        self.n_lines = len(self.lines)
        self.n_gens = len(self.generators)

        self.b = None
        self.g = None
        self.G = None
        self.B = None

        self.load_flow_solved = False
        self.Ybus_load = None
    
    @property
    def N(self):
        return len(self.buses)

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
    def from_excel(filename):
        pass

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
                self.n_buses = len(self.buses)
            elif object_class is Line:
                self.lines = collection
                self.n_lines = len(self.lines)
            elif object_class is Generator:
                self.generators = collection
                self.n_gens = len(self.generators)
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
        Ya = self.Ybus[:self.n_gens, :self.n_gens]

        # Yb is a matrix of size (n_gens, n_buses-n_gens)
        Yb = self.Ybus[:self.n_gens, self.n_gens:]

        # Yc is the transpose of Yb
        Yc = Yb.T

        # Yd is a matrix that goes from n_gens: and n_gens:
        Yd = self.Ybus[self.n_gens:, self.n_gens:]

        # Now, calculate Kron eduction
        Ykron = Ya - Yb@(np.linalg.inv(Yd)@Yc)

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
    
