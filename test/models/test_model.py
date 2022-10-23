# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Tests of the actual solver model and results.
"""
import os
import copy
import unittest
import numpy as np

from lindbladmpo.LindbladMPOSolver import LindbladMPOSolver
from lindbladmpo.examples.simulation_building.LindbladMatrixSolver import (
    LindbladMatrixSolver,
)

s_output_path = (
    os.path.abspath("./output") + "/"
)  # All solver output files will be written here
if not os.path.exists(s_output_path):
    os.mkdir(s_output_path)
s_cygwin_path = None
s_solver_path = None


class LindbladMPOSolverModel(unittest.TestCase):
    """The class testing the actual solver model and numerical results."""

    def test_all_zero(self):
        """Test a trivial setup."""
        solver_params = {
            "tau": 0.1,
            "t_final": 1,
            "N": 2,
            "g_1": 0,
            "output_files_prefix": s_output_path + "test_all_zero",
            "1q_components": ["X", "Y", "Z"],
            "b_quiet": True,
        }
        solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
        solver.solve()
        expected_XY = 0
        expected_Z = 1
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (0,))][1][-1], expected_Z)
        # second qubit should be the same
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (1,))][1][-1], expected_Z)

    def test_hz_nonzero(self):
        """Test the dynamics with h_z nonzero."""
        solver_params = {
            "tau": 0.1,
            "t_final": 1,
            "N": 2,
            "g_1": 0,
            "l_x": 0,
            "h_z": 1,
            "output_files_prefix": s_output_path + "test_hz_nonzero",
            "1q_components": ["X", "Y", "Z"],
            "b_quiet": True,
        }
        solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
        solver.solve()
        expected_XY = 0
        expected_Z = 1
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (0,))][1][-1], expected_Z)
        # second qubit should be the same
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (1,))][1][-1], expected_Z)

    def test_steady_state(self):
        """Test a strong decay to the Bloch sphere's south pole."""
        solver_params = {
            "tau": 0.1,
            "t_final": 10,
            "N": 2,
            "g_1": 2,
            "l_x": 0,
            "h_z": 1,
            "output_files_prefix": s_output_path + "test_steady_state",
            "1q_components": ["X", "Y", "Z"],
            "b_quiet": True,
        }
        solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
        solver.solve()
        expected_XY = 0
        expected_Z = -1
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (0,))][1][-1], expected_Z)
        # second qubit should be the same
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (1,))][1][-1], expected_Z)

    def test_steady_state_2(self):
        """Test an intermediate-z steady state."""
        solver_params = {
            "tau": 0.1,
            "t_final": 10,
            "N": 2,
            "g_1": 5,
            "g_0": 1,
            "l_x": 0,
            "h_z": 5,
            "output_files_prefix": s_output_path + "test_steady_state_2",
            "1q_components": ["X", "Y", "Z"],
            "b_quiet": True,
        }
        solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
        solver.solve()
        expected_XY = 0
        expected_Z = -4 / 6
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (0,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (0,))][1][-1], expected_Z)
        # second qubit should be the same
        self.assertAlmostEqual(solver.result["obs-1q"][("x", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("y", (1,))][1][-1], expected_XY)
        self.assertAlmostEqual(solver.result["obs-1q"][("z", (1,))][1][-1], expected_Z)

    def test_3Q_two_solvers(self):
        """Test a trivial setup."""
        N = 3
        t_final = 3
        J = np.zeros((N, N))
        J[0, 1] = -0.5
        J[2, 1] = 0.8
        J_z = np.zeros((N, N))
        J_z[1, 2] = 0.7
        J_z[0, 2] = -0.4
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "g_0": [0.1, 0.5, 0.2],
            "g_1": 0.2,
            "g_2": 0.1,
            "h_x": [0.6, 1.2, 0.7],
            "h_y": [-0.6, 0.2, 0.9],
            "h_z": [0.1, -1.2, 0.7],
            "J": J,
            "J_z": J_z,
            "init_product_state": ["+x", 0.8, (0.3, -1.8)],
            "init_cz_gates": [(0, 2)],
            "1q_components": ["X", "Y", "Z"],
            "2q_components": ["XY", "XZ", "ZZ"],
            "b_save_final_state": True,
            "b_quiet": True,
        }
        solver_params2 = copy.deepcopy(solver_params1)

        s_files_prefix1 = s_output_path + "test_3Q_two_solvers_1"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        s_files_prefix2 = s_output_path + "test_3Q_two_solvers_2"
        solver_params2.update({"output_files_prefix": s_files_prefix2})
        solver2 = LindbladMatrixSolver(solver_params1)
        solver2.solve()
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (0,))][1][-1],
            solver2.result["obs-1q"][("x", (0,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (1, 2))][1][-1],
            solver2.result["obs-2q"][("xz", (1, 2))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("zz", (2, 0))][1][-1],
            solver2.result["obs-2q"][("zz", (2, 0))][1][-1],
        )

        load_dict = {
            "b_save_final_state": False,
            "init_product_state": None,
            "init_cz_gates": [],
            "t_init": t_final,
            "t_final": t_final + 2,
        }
        solver_params1.update(load_dict)
        solver_params1.update(
            {
                "load_files_prefix": s_files_prefix1,
                "output_files_prefix": s_files_prefix1 + "B",
            }
        )
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()
        solver_params2.update(load_dict)
        solver_params2.update(
            {
                "load_files_prefix": s_files_prefix2,
                "output_files_prefix": s_files_prefix2 + "B",
            }
        )
        solver2 = LindbladMatrixSolver(solver_params1)
        solver2.solve()
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("y", (1,))][1][-1],
            solver2.result["obs-1q"][("y", (1,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("z", (2,))][1][-1],
            solver2.result["obs-1q"][("z", (2,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (0, 1))][1][-1],
            solver2.result["obs-2q"][("xy", (0, 1))][1][-1],
        )


if __name__ == "__main__":
    unittest.main()
