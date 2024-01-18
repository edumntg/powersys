from src.powersys import *
from src.powersys.models import PowerSystem, PowerSystemArgs
from src.powersys.solvers import LFSolver
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

import unittest

class TestKronReduction(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        #  Load system
        args = PowerSystemArgs(
            f = 60,
            buses = PowerSystem.load_buses(CWD + '/sample_data/ieee9_buses.csv'),
            lines = PowerSystem.load_lines(CWD + '/sample_data/ieee9_lines.csv'),
            generators = PowerSystem.load_gens(CWD + '/sample_data/ieee9_gens.csv')
        )

        system = PowerSystem(args)

        solver = LFSolver(system)
        solver.solve(disp = False)

        # Now, construct the Ybus-load
        system.construct_load_ybus()

        # Construct Kron
        system.kron_reduction()

        # Construct Yrm
        system.YRM()

        # Construct RM vectors
        system.rm()

        system.compute_terminal_values()

        system.compute_gap_power()

        self.system = system

    def test_correct_value(self):
        expected = np.array([
            [1.1051-4.6957*1j, 0.0965+2.2570*1j, 0.0046+2.2748*1j],
            [0.0965+2.2570*1j, 0.7355-5.1143*1j, 0.1230+2.8257*1j],
            [0.0046+2.2748*1j, 0.1230+2.8257*1j, 0.7214-5.0231*1j]
        ])

        actual = self.system.Ykron

        self.assertTrue(np.abs(np.subtract(expected, actual)).max() <= 1E-3)

    def test_correct_shape(self):
        self.assertTrue(self.system.Ykron.shape == (self.system.n_gens, self.system.n_gens))

if __name__ == '__main__':
    unittest.main()
