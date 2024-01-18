from src.powersys import *
from src.powersys.solvers import LFSolver
from src.powersys.models import PowerSystem, PowerSystemArgs
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

import unittest

class TestGapPowerCalculation(unittest.TestCase):

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
        expected = np.array([0.7170, 1.6325, 0.8521])
        actual = np.array([gen.Pm for gen in self.system.generators])

        self.assertTrue(np.abs(np.subtract(expected, actual)).max() <= 1E-3)

    def test_correct_size(self):
        self.assertTrue(np.array([gen.Pm for gen in self.system.generators]).size == self.system.n_gens)

if __name__ == '__main__':
    unittest.main()

