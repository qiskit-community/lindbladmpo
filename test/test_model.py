# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from MPO_Lindblad_solver.MPOLindbladSolver import MPOLindbladSolver
import unittest
import numpy as np

s_path = "C:/Users/galvz/PycharmProjects/sim_func/"
s_executable = "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator/lindblad.exe"


class TestModel(unittest.TestCase):
    def test_All_zero(self):
        solver_params = {'tau': 1, 't_final': 1, 'N': 1, 'g_1': 0,
                         'input_file_prefix': s_path + "MPO.input",
                         'output_file_prefix': s_path + "MPO"}
        solver = MPOLindbladSolver(solver_params, "C:/cygwin64/bin/bash.exe", s_executable)
        solver.solve()
        expected_XY = 0
        expected_Z = 1
        self.assertEqual(solver.result['1q'][(1, 'X', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Y', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Z', 1)], expected_Z)


    def test_hz_not_zero(self):
        solver_params = {'tau': 1, 't_final': 1, 'N': 1, 'g_1': 0, 'l_x': 0, 'h_z': 5,
                         'input_file_prefix': s_path + "MPO.input",
                         'output_file_prefix': s_path + "MPO"}
        solver = MPOLindbladSolver(solver_params, "C:/cygwin64/bin/bash.exe", s_executable)
        solver.solve()
        expected_XY = 0
        expected_Z = 1
        self.assertEqual(solver.result['1q'][(1, 'X', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Y', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Z', 1)], expected_Z)

    def test_steady_state(self):
        solver_params = {'tau': 1, 't_final': 1, 'N': 1, 'g_1': 5, 'l_x': 0, 'h_z': 5,
                         'input_file_prefix': s_path + "MPO.input",
                         'output_file_prefix': s_path + "MPO"}
        solver = MPOLindbladSolver(solver_params, "C:/cygwin64/bin/bash.exe", s_executable)
        solver.solve()
        expected_XY = 0
        expected_Z = 1
        self.assertEqual(solver.result['1q'][(1, 'X', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Y', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Z', 1)], expected_Z)

    def test_steady_state_2(self):
        solver_params = {'tau': 1, 't_final': 10, 'N': 1, 'g_1': 5, 'g_0': 1, 'l_x': 0, 'h_z': 5,
                         'input_file_prefix': s_path + "MPO.input",
                         'output_file_prefix': s_path + "MPO"}
        solver = MPOLindbladSolver(solver_params, "C:/cygwin64/bin/bash.exe", s_executable)
        solver.solve()
        expected_XY = 0
        expected_Z = 4/6
        self.assertEqual(solver.result['1q'][(1, 'X', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Y', 1)], expected_XY)
        self.assertEqual(solver.result['1q'][(1, 'Z', 1)], expected_Z)


if __name__ == '__main__':
    unittest.main()
