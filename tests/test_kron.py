import powersys as ps
import pandas as pd
import numpy as np
import os

CWD = os.path.dirname(os.path.abspath(__file__))

import unittest

class TestKronReduction(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        #  Load system
        args = ps.models.PowerSystemArgs(
            f = 60,
            buses = ps.models.PowerSystem.load_buses(CWD + '/../sample_data/ieee9_buses.csv'),
            lines = ps.models.PowerSystem.load_lines(CWD + '/../sample_data/ieee9_lines.csv'),
            generators = ps.models.PowerSystem.load_gens(CWD + '/../sample_data/ieee9_gens.csv')
        )

        system = ps.models.PowerSystem(buses = args.buses, lines = args.lines, generators = args.generators)

        solver = ps.solvers.LF(system)
        solver.solve(disp = False, method = "newton-raphson")

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
