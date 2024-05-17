from ..objects import *

class DataLoader(object):
    
    def __init__(self, filename = '', buses = ObjectCollection(), lines = ObjectCollection(), generators = ObjectCollection()):
        self.filename = filename
        self.buses = buses
        self.lines = lines
        self.generators = generators

    def load(self):
        assert not self.filename is None, "Invalid filename"

        extension = self.filename.split(".")[1]
        if extension == "txt":
            new_loader = DataLoader.from_text(self.filename)
            self.buses = new_loader.buses
            self.lines = new_loader.lines
            self.generators = new_loader.generators

    def __getitem__(self, idx):
        return self.buses[idx], self.lines[idx], self.generators[idx]

    @staticmethod
    def from_text(filename):
        """ NOTE: In order for this function to work, the .txt file must be in IEEE format """

        buses = ObjectCollection()
        ps_lines = ObjectCollection()
        generators = ObjectCollection()

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

            buses.add(Busbar(
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

            ps_lines.add(Line(
                from_bus = from_bus,
                to_bus = to_bus,
                R = R,
                X = X,
                B = B
            ))

        return DataLoader(filename, buses, ps_lines, generators)
