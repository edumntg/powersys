import powersys as ps
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

import unittest

class TestJSONLoader(unittest.TestCase):

    def test_read_json(self):
        system = ps.models.PowerSystem.read_json(CWD + '/sample_data/ieee9.json')

        assert system.N == 9, "Incorrect number of bus bars"
        assert system.n_gens == 3, "Incorrect number of generators"
        assert system.n_lines == 9, "Incorrect number of lines"
        assert system.f == 60, "Incorrect frequency value"

if __name__ == '__main__':
    unittest.main()