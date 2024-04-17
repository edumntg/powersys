import pandas as pd
from . import Bus

class ObjectCollection(object):

    def __init__(self):
        self.items = []

    def add(self, item):
        self.items.append(item)

    def append(self, item):
        return self.add(item)

    def __getitem__(self, index):
        if index < 0 or index >= len(self.items):
            raise IndexError
        
        return self.items[index]
    
    def __str__(self):
        return "\n".join([str(item) for item in self.items])
    
    def __len__(self):
        return len(self.items)
    
    @staticmethod
    def from_csv(filename, object_class, delimiter = ','):
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
    
    def from_dict(self, dictio, object_class = None):
        assert object_class is not None, "object_class property cannot be null"
        assert isinstance(dictio, dict), "Data is not a valid dictionary"

        self.items = []

        for (key, value) in dictio.items():
            object_id = int(key)
            dict_data = dict(value)
            dict_data['id'] = object_id

            self.items.append(object_class.from_dict(dict_data))

        return self



