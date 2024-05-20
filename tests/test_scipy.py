import powersys as ps
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

import unittest

class TestLFScipy(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        #  Load system
        args = ps.models.PowerSystemArgs(
            f = 60,
            buses = ps.models.PowerSystem.load_buses(CWD + '/../sample_data/ieee9_buses.csv'),
            lines = ps.models.PowerSystem.load_lines(CWD + '/../sample_data/ieee9_lines.csv'),
            generators = ps.models.PowerSystem.load_gens(CWD + '/../sample_data/ieee9_gens.csv')
        )

        self.system = ps.models.PowerSystem(buses = args.buses, lines = args.lines, generators = args.generators)

    def test_lf(self):
        solver = ps.solvers.LF(self.system)
        try:
            solver.solve(disp = False, method = "fsolve")
            self.assertTrue(True)
        except:
            self.assertFalse(True)
            
if __name__ == '__main__':
    unittest.main()