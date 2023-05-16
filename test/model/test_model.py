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
from lindbladmpo.plot_routines import (
    prepare_2q_density_operator,
    prepare_concurrence_data,
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
        """Test trivial dynamics with h_z nonzero."""
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
        """Test a steady state with an intermediate z value."""
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

    def test_product_state(self):
        """Test multiple product state initializations that are equivalent."""
        N = 6
        t_final = 1.5
        J = np.zeros((N, N))
        J[0, 1] = -1.8
        J[2, 3] = -1.8
        J[4, 5] = -1.8
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "g_1": 0.1,
            "J": J,
            "init_product_state": [
                "+y",
                "+x",
                "+y",
                (np.pi / 2, 0.0),
                (0.5, 0.0, -0.5),
                (0.5, 0.5, 0.0),
            ],
            "1q_components": ["X", "Y"],
            "2q_components": ["XY"],
            "b_quiet": True,
        }
        s_files_prefix1 = s_output_path + "test_product_state"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (1,))][1][0],
            solver1.result["obs-1q"][("x", (3,))][1][0],
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (3,))][1][0],
            solver1.result["obs-1q"][("x", (5,))][1][0],
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (3,))][1][0],
            1.0,
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("y", (4,))][1][0],
            1.0,
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (1,))][1][-1],
            solver1.result["obs-1q"][("x", (3,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("y", (2,))][1][-1],
            solver1.result["obs-1q"][("y", (4,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xy", (2, 3))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (2, 3))][1][-1],
            solver1.result["obs-2q"][("xy", (4, 5))][1][-1],
        )

    def test_cz_gates(self):
        """Test multiple product state initializations that are equivalent."""
        N = 6
        t_final = 1.5
        J_z = np.zeros((N, N))
        J_z[0, 1] = 1.8
        J_z[2, 3] = 1.8
        J_z[4, 5] = 1.8
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "g_1": 0.1,
            "J_z": J_z,
            "init_product_state": [
                "id",
                "+x",
                0.5,
                (np.pi / 2, 0.0),
                (0.5, 0.0, 0.0),
                (0.5, 0.5, 0.0),
            ],
            "init_cz_gates": [(0, 1), (2, 3), (4, 5)],
            "1q_components": ["X", "Y"],
            "2q_components": ["XZ"],
            "b_quiet": True,
        }
        s_files_prefix1 = s_output_path + "test_cz_gates"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (1, 0))][1][0],
            1.0,
        )
        self.assertAlmostEqual(
            solver1.result["global"][("osee_center", ())][1][0],
            np.log(2),
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (0,))][1][-1],
            solver1.result["obs-1q"][("x", (2,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("y", (3,))][1][-1],
            solver1.result["obs-1q"][("y", (5,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (1, 0))][1][-1],
            solver1.result["obs-2q"][("xz", (3, 2))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (3, 2))][1][-1],
            solver1.result["obs-2q"][("xz", (5, 4))][1][-1],
        )

    def test_3Q_two_solvers(self):
        """Compare the very generic evolution of three qubits using two solvers.
        The two solvers employed are the MPO solver and the scipy solver using qiskit-dynamics.
        Moreover, the evolution is saved at some intermediate time and loaded as the initial
        condition for some further time propagation, and then compared again.
        """
        N = 3
        t_final = 3
        t_delta = 2
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
            "init_product_state": ["+x", (0.8, 0.1, -0.6), (0.3, -1.8)],
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
        solver2 = LindbladMatrixSolver(solver_params2)
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
            "t_final": t_final + t_delta,
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
        solver2 = LindbladMatrixSolver(solver_params2)
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

    def test_3Q_again(self):
        """Compare the very generic evolution of three qubits using two solvers."""
        N = 3
        t_final = 2
        t_delta = 2
        J = np.zeros((N, N))
        J[0, 1] = 1.5
        J[2, 1] = 0.8
        J_z = np.zeros((N, N))
        J_z[1, 2] = 0.7
        J_z[0, 2] = 1.4
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "g_0": [0.1, 0.05, 0.02],
            "g_1": 0.02,
            "g_2": 0.1,
            "h_x": [0.6, 1.2, 0.7],
            "h_y": [-0.6, 0.2, 0.9],
            "h_z": [0.1, -1.2, 0.7],
            "J": J,
            "J_z": J_z,
            "init_product_state": ["id", 0.2, (0.2, 0.0, -0.9)],
            "init_cz_gates": [(1, 2)],
            "1q_components": ["X", "Y", "Z"],
            "2q_components": ["XY", "YX", "XX"],
            "3q_components": ["XYZ", "ZYX", "XZX"],
            "3q_indices": [(0, 1, 2)],
            "b_save_final_state": True,
            "b_quiet": True,
        }
        solver_params2 = copy.deepcopy(solver_params1)

        s_files_prefix1 = s_output_path + "test_3Q_again_1"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        s_files_prefix2 = s_output_path + "test_3Q_again_2"
        solver_params2.update({"output_files_prefix": s_files_prefix2})
        solver2 = LindbladMatrixSolver(solver_params2)
        solver2.solve()
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (0,))][1][-1],
            solver2.result["obs-1q"][("x", (0,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xx", (1, 2))][1][-1],
            solver2.result["obs-2q"][("xx", (2, 1))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (2, 0))][1][-1],
            solver2.result["obs-2q"][("yx", (0, 2))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-3q"][("xyz", (0, 1, 2))][1][-1],
            solver2.result["obs-3q"][("xyz", (0, 1, 2))][1][-1],
            places=5,
        )
        self.assertAlmostEqual(
            solver1.result["obs-3q"][("zyx", (0, 1, 2))][1][-1],
            solver2.result["obs-3q"][("zyx", (0, 1, 2))][1][-1],
            places=5,
        )
        self.assertAlmostEqual(
            solver1.result["obs-3q"][("xzx", (0, 1, 2))][1][-1],
            solver2.result["obs-3q"][("xzx", (0, 1, 2))][1][-1],
            places=5,
        )
        solver_params1A = solver_params1.copy()
        solver_params1A.update(
            {
                "t_final": t_final + t_delta,
                "output_files_prefix": s_files_prefix1 + "A",
            }
        )
        solver1A = LindbladMPOSolver(solver_params1A, s_cygwin_path, s_solver_path)
        solver1A.solve()
        solver_params2A = solver_params2.copy()
        solver_params2A.update(
            {
                "t_final": t_final + t_delta,
                "output_files_prefix": s_files_prefix2 + "A",
            }
        )
        solver2A = LindbladMPOSolver(solver_params2A, s_cygwin_path, s_solver_path)
        solver2A.solve()

        load_dict = {
            "b_save_final_state": False,
            "init_product_state": None,
            "init_cz_gates": [],
            "t_init": t_final,
            "t_final": t_final + t_delta,
        }
        solver_params1.update(load_dict)
        solver_params1.update(
            {
                "load_files_prefix": s_files_prefix1,
                "output_files_prefix": s_files_prefix1 + "B",
            }
        )
        solver1B = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1B.solve()
        solver_params2.update(load_dict)
        solver_params2.update(
            {
                "load_files_prefix": s_files_prefix2,
                "output_files_prefix": s_files_prefix2 + "B",
            }
        )
        solver2B = LindbladMatrixSolver(solver_params2)
        solver2B.solve()
        self.assertAlmostEqual(
            solver1A.result["obs-1q"][("y", (1,))][1][-1],
            solver2B.result["obs-1q"][("y", (1,))][1][-1],
            6,
        )
        self.assertAlmostEqual(
            solver1A.result["obs-1q"][("z", (2,))][1][-1],
            solver2B.result["obs-1q"][("z", (2,))][1][-1],
            6,
        )
        self.assertAlmostEqual(
            solver1A.result["obs-2q"][("xy", (0, 1))][1][-1],
            solver2B.result["obs-2q"][("xy", (0, 1))][1][-1],
            6,
        )

    def test_dm(self):
        """Test density matrix reconstruction."""
        N = 2
        t_final = 1
        solver_params1 = {
            "tau": 0.1,
            "t_final": t_final,
            "N": N,
            "init_product_state": [
                "+x",
                "id",
            ],
            "init_cz_gates": [(0, 1)],
            "1q_components": ["x", "y", "z"],
            "2q_components": ["xx", "yy", "zz", "xy", "xz", "yz"],
            "b_quiet": True,
        }
        s_files_prefix1 = s_output_path + "test_dm"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()
        (_, rho_list), _ = prepare_2q_density_operator(solver1.result, (0, 1), [-1])

        expected_rho = 0.25 * np.asarray(
            np.diag([1, 1, 1, 1]), dtype=complex
        ) + 0.25 * np.kron(
            np.asarray([[0, 1], [1, 0]], dtype=complex),
            np.asarray([[1, 0], [0, -1]], dtype=complex),
        )
        self.assertAlmostEqual(
            np.linalg.norm(rho_list[0] - expected_rho),
            0.0,
        )

    def test_graph_state_concurrence(self):
        """Test concurrence calculation."""
        N = 2
        t_final = 1
        solver_params1 = {
            "tau": 0.1,
            "t_final": t_final,
            "N": N,
            "init_graph_state": [(0, 1)],
            "1q_components": ["x", "y", "z"],
            "2q_components": ["xx", "yy", "zz", "xy", "xz", "yz"],
            "b_quiet": True,
        }
        s_files_prefix1 = s_output_path + "test_graph_state_concurrence"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()
        (_, c_data), _ = prepare_concurrence_data(solver1.result, (0, 1))
        self.assertAlmostEqual(
            c_data[0],
            1.0,
        )

    def test_apply_gates_1(self):
        """Compare the application of gates at intermediate times."""
        N = 2
        t_final = 2
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "h_z": 0.38,
            "J_z": 2 * np.pi / t_final,
            "1q_components": ["X", "Y"],
            "2q_components": ["XZ", "ZZ"],
            "b_quiet": True,
            "init_product_state": ["-x", "+y"],
            "apply_gates": [(t_final / 2.0, "x", 0), (t_final / 2.0, "y", 1)],
        }
        s_files_prefix1 = s_output_path + "test_apply_gates_1A"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        # The qubits should disentangle (since the time is exactly a period of J_z),
        # and return to their initial state since the x and y gates are applied at half
        # the simulation time, doing a hahn echo that cancels the h_z rotation
        self.assertAlmostEqual(solver1.result["obs-1q"][("x", (0,))][1][-1], -1.0)
        self.assertAlmostEqual(solver1.result["obs-1q"][("y", (1,))][1][-1], 1.0)
        self.assertAlmostEqual(solver1.result["obs-2q"][("xz", (0, 1))][1][-1], 0.0)
        self.assertAlmostEqual(solver1.result["obs-2q"][("xz", (1, 0))][1][-1], 0.0)

    def test_apply_gates_2(self):
        """Compare the application of gates at intermediate times."""
        N = 2
        t_final = 1
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "g_2": 0.1,
            "h_z": -0.2,
            "J_z": -0.8,
            "1q_components": ["X", "Y"],
            "2q_components": ["XZ", "ZZ"],
            "b_quiet": True,
        }
        solver_params2 = copy.deepcopy(solver_params1)
        solver_params1.update(
            {
                "init_product_state": ["-x", "+y"],
                "init_cz_gates": [(0, 1)],
                "apply_gates": [(0.6001, "x", 1), (0.599, "cz", 1, 0)],
            }
        )
        solver_params2.update(
            {
                "init_product_state": ["+x", "-y"],
                "apply_gates": [
                    (0.0, "y", 0),
                    (0.0, "z", 1),
                    (0.0, "cz", 0, 1),
                    (0.60, "x", 1),
                    (0.60, "cz", 1, 0),
                ],
            }
        )
        s_files_prefix1 = s_output_path + "test_apply_gates_2A"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        s_files_prefix2 = s_output_path + "test_apply_gates_2B"
        solver_params2.update({"output_files_prefix": s_files_prefix2})
        solver2 = LindbladMPOSolver(solver_params2, s_cygwin_path, s_solver_path)
        solver2.solve()
        self.assertAlmostEqual(
            solver1.result["obs-1q"][("x", (0,))][1][-1],
            solver2.result["obs-1q"][("x", (0,))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (0, 1))][1][-1],
            solver2.result["obs-2q"][("xz", (0, 1))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (1, 0))][1][-1],
            solver2.result["obs-2q"][("xz", (1, 0))][1][-1],
        )

    def test_apply_gates_cx_sx_h(self):
        """Compare the application of gates at intermediate times."""
        N = 8
        t_final = 1
        solver_params1 = {
            "tau": 0.02,
            "t_final": t_final,
            "N": N,
            "h_z": 0.43,
            "1q_components": ["x", "y", "z"],
            "2q_components": ["xx", "yy", "zz", "xy", "xz", "yz"],
            "b_quiet": True,
            "init_product_state": ["+x", "+x", "+z", "+z", "+z", "+z", "-z", "-z"],
            "apply_gates": [
                (0.0, "cz", 0, 1),
                (0.0, "h", 2),
                (0.0, "h", 3),
                (0.0, "cz", 2, 3),
                (0.0, "h", 4),
                (0.0, "h", 5),
                (0.0, "h", 5),
                (0.0, "cx", 4, 5),
                (0.0, "h", 5),
                (0.0, "sx", 6),
                (0.0, "sx", 7),
                (0.0, "cz", 6, 7),
            ],
        }
        s_files_prefix1 = s_output_path + "test_apply_gates_cx_sx_h"
        solver_params1.update({"output_files_prefix": s_files_prefix1})
        solver1 = LindbladMPOSolver(solver_params1, s_cygwin_path, s_solver_path)
        solver1.solve()

        # Below we test exact equality of 2Q observables of the first three pairs -
        # which should all have identical states, and we test that the fourth pair
        # has a concurrence equal to the first pair's (would be 1), since it's not
        # identical to it (being prepared with an sx gate rather than h).
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xz", (2, 3))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xx", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xx", (2, 3))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xy", (2, 3))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xz", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xz", (4, 5))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xx", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xx", (4, 5))][1][-1],
        )
        self.assertAlmostEqual(
            solver1.result["obs-2q"][("xy", (0, 1))][1][-1],
            solver1.result["obs-2q"][("xy", (4, 5))][1][-1],
        )

        (_, c_data1), _ = prepare_concurrence_data(solver1.result, (0, 1))
        (_, c_data2), _ = prepare_concurrence_data(solver1.result, (6, 7))
        self.assertAlmostEqual(c_data1[-1], 1.0)
        self.assertAlmostEqual(c_data1[-1], c_data2[-1])


if __name__ == "__main__":
    unittest.main()
