# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import math
import time
from typing import Sized, Iterable, Container
import scipy

from lindbladmpo.LindbladMPOSolver import *
from .operators import *
from .operators_library import *

from qiskit.quantum_info import DensityMatrix
from qiskit_dynamics.array import Array
from qiskit_dynamics import solve_lmde
from qiskit_dynamics.signals import Signal
from qiskit_dynamics.models import HamiltonianModel, LindbladModel


class LindbladMatrixSolver(LindbladMPOSolver):
    """A solver duplicating the interface of LindbladMPOSolver, but uses qiskit-dynamics.

    This class makes it possible to define a simulation problem that can be solved either using
    the MPO solver or using qiskit-dynamics with a very similar interface. It generates identical
    input and output files, making it possible to use all the analysis and plotting functions
    in the lindbladmpo package in a transparent way.
    """

    IMAGINARY_THRESHOLD = 1e-4
    """Threshold for the imaginary value of a quantity that should be real, to issue a warning."""

    TRACE_RHO_WARN_THRESHOLD = 1e-4
    """Threshold for the deviation of the density matrix trace from 1, to issue a warning."""

    DEFAULT_PARAMETERS = {
        "t_init": (0.0, "s"),
        "output_files_prefix": ("lindblad", "s"),
        "b_unique_id": (False, "s"),
        "h_x": (0.0, "v"),
        "h_y": (0.0, "v"),
        "h_z": (0.0, "v"),
        "J_z": (0.0, "m"),
        "J": (0.0, "m"),
        "g_0": (0.0, "v"),
        "g_1": (0.0, "v"),
        "g_2": (0.0, "v"),
        "g_3": (0.0, "v"),
        "g_4": (0.0, "v"),
        "b_quiet": (False, "s"),
        "b_save_final_state": (False, "s"),
        "1q_components": (["z"], "s"),
        "2q_components": (["zz"], "s"),
        "method": ("RK45", "s"),
        "atol": (1e-8, "s"),
        "rtol": (1e-8, "s"),
    }
    """Default values for some parameters, with an indication of being a scalar, vector, or matrix."""

    def __init__(self, parameters: Optional[dict] = None):
        """Initialize the instance, and possibly build the input file from the parameters.

        Args:
                parameters: The model parameters. If not None, the member build() function is
                        invoked to build the input file.
        """
        super().__init__(parameters)
        self._log_file = None

    def solve(self):
        """Solves the simulation and loads the result dictionaries."""
        if self.s_input_file == "":
            self.build()

        file_1q = None
        file_2q = None
        file_3q = None
        file_cu = None
        file_gl = None
        self._log_file = open(self.s_output_path + ".log.txt", "w")
        try:
            input_file = open(self.s_input_file, "r")
            self._print(
                "-------------------------------------------------------------------------"
            )
            self._print(input_file.read())
            self._print(
                "-------------------------------------------------------------------------"
            )
            self._print(
                "Creating solver matrices and executing scipy solver using qiskit-dynamics."
            )
            ts0 = time.time()
            parameters = self.parameters
            n_qubits = parameters.get("N", None)
            t_final = parameters.get("t_final", None)
            tau = parameters.get("tau", None)
            if n_qubits is None or t_final is None or tau is None:
                raise Exception(
                    "The three input parameters 'N', 't_final', and 'tau' must be assigned, "
                    "as they do not have a default value."
                )
            s_load_files_prefix: str = parameters.get("load_files_prefix", "")

            init_graph_state = self._get_parameter("init_graph_state")
            init_cz_gates = self._get_parameter("init_cz_gates")
            b_graph_state, b_cz_pairs = False, False
            s_cz_param = ""
            if init_graph_state is not None and len(init_graph_state) > 0:
                b_cz_pairs = True
                b_graph_state = True
                s_cz_param = "init_graph_state"
                if init_cz_gates is None or len(init_cz_gates) == 0:
                    init_cz_gates = init_graph_state
                else:
                    raise Exception(
                        "The parameter init_cz_gates cannot be used "
                        "if init_graph_state is nonempty."
                    )
            elif init_cz_gates is not None and len(init_cz_gates) > 0:
                b_cz_pairs = True
                s_cz_param = "init_cz_gates"

            s_init_param = ""
            init_pauli_state = self._get_parameter("init_pauli_state")
            init_product_state: Any = self._get_parameter("init_product_state")
            if init_pauli_state is not None and (
                self.is_float(init_pauli_state) or len(init_pauli_state) != 0
            ):
                s_init_param = "init_pauli_state"
                self._print(
                    "Warning: the init_pauli_state parameter has been deprecated and "
                    "will be removed in the future. Please use init_product_state instead."
                )
                if init_product_state is None or (
                    not self.is_float(init_product_state)
                    and len(init_product_state) == 0
                ):
                    init_product_state = init_pauli_state
                else:
                    raise Exception(
                        "The parameter init_pauli_state cannot be used "
                        "if init_product_state is nonempty."
                    )
            else:
                s_init_param = "init_product_state"
            if init_product_state is None:
                init_len = 0
            elif isinstance(init_product_state, str) and len(init_product_state) == 0:
                init_len = 0
            else:
                if (
                    self.is_float(init_product_state)
                    or isinstance(init_product_state, str)
                    or isinstance(init_product_state, tuple)
                ):
                    init_product_state = [init_product_state]
                init_len = len(init_product_state)

            if s_load_files_prefix != "" and (b_cz_pairs or init_len > 0):
                raise Exception(
                    "If load_files_prefix is nonempty, no other initialization "
                    "parameter can be used"
                )
            if init_len == 0:
                if b_graph_state:
                    init_product_state = ["+x"] * n_qubits
                else:
                    init_product_state = ["+z"] * n_qubits
            else:
                if b_graph_state:
                    raise Exception(
                        "If init_graph_state is nonempty, no other initialization "
                        "parameter can be used."
                    )
                if init_len == 1:
                    init_product_state = [init_product_state[0]] * n_qubits
                elif init_len != n_qubits:
                    raise Exception(
                        f"The parameter {s_init_param} has {init_len}"
                        f"value(s) but 0, 1 or {n_qubits} value(s) were expected."
                    )

            in_cz_gate = np.zeros((n_qubits,))
            if init_cz_gates is not None:
                for q_pair in init_cz_gates:
                    i = q_pair[0]
                    j = q_pair[1]
                    if i >= n_qubits or j >= n_qubits:
                        raise Exception(
                            f"The parameter {s_cz_param} contains an invalid index pair, ({i}, {j})."
                        )
                    in_cz_gate[i] = 1
                    in_cz_gate[j] = 1

            apply_gates = self._get_parameter("apply_gates")
            if apply_gates is not None and len(apply_gates):
                raise Exception(
                    "The parameter apply_gates is not currently supported by "
                    "this solver. Please use the MPO solver for this simulation."
                )

            h_x = self._get_parameter("h_x")
            h_y = self._get_parameter("h_y")
            h_z = self._get_parameter("h_z")
            g_0 = self._get_parameter("g_0")
            g_1 = self._get_parameter("g_1")
            g_2 = self._get_parameter("g_2")
            g_3 = self._get_parameter("g_3")
            g_4 = self._get_parameter("g_4")
            J = self._get_parameter("J")
            J_z = self._get_parameter("J_z")
            _1q_components = self._get_parameter("1q_components")
            _1q_indices = parameters.get("1q_indices", None)
            if _1q_indices is None:  # Add 1Q observables for all qubits.
                _1q_indices = range(0, n_qubits)

            _2q_components = self._get_parameter("2q_components")
            _2q_indices = parameters.get("2q_indices", None)
            if _2q_indices is None:  # Add 2Q observables for all qubit pairs.
                _2q_indices = []
                for i in range(0, n_qubits):
                    for j in range(0, n_qubits):
                        if i != j:
                            _2q_indices.append((i, j))

            _3q_components = parameters.get("3q_components", None)
            if _3q_components is None:  # No observables are added
                _3q_components = []
            _3q_indices = parameters.get("3q_indices", None)
            if _3q_indices is None:  # No observables are added
                _3q_indices = []
            custom_observables = parameters.get("custom_observables", None)
            if custom_observables is None:  # No observables are added
                custom_observables = []

            t_init = self._get_parameter("t_init")
            method = self._get_parameter("method")
            atol = self._get_parameter("atol")
            rtol = self._get_parameter("rtol")

            r_qubits = range(n_qubits)
            subsystem_dims = OrderedDict()
            H = 0.0 * Id(0)  # just dummy initialization
            rho_0 = Id(0)  # just dummy initialization
            gs = Id(0)  # just dummy initialization
            id_op = Id(0)  # just dummy initialization

            L_ops = []
            L_sig = []
            obs_1q = []
            obs_2q = []
            obs_3q = []
            obs_1q_key = []
            obs_2q_key = []
            obs_3q_key = []
            obs_cu_mat = []
            obs_cu_key = []
            for i_qubit in r_qubits:
                subsystem_dims[i_qubit] = 2
                gs *= PlusZ(i_qubit)
                id_op *= Id(i_qubit)
                if s_load_files_prefix == "":
                    q_init: Any = init_product_state[i_qubit]
                    b_tuple = isinstance(q_init, tuple)
                    b_diagonal = self.is_float(q_init) or (b_tuple and len(q_init) == 1)
                    if b_diagonal:
                        b = q_init[0] if b_tuple else q_init
                        diagonal = [b, 1.0 - b]
                        rho_0 *= Diagonal(i_qubit, diagonal)
                    elif b_tuple:
                        if len(q_init) == 2:
                            theta: float = q_init[0]
                            phi: float = q_init[1]
                            rho_0 *= PolarState(i_qubit, theta, phi)
                        elif len(q_init) == 3:
                            a: float = q_init[0]
                            b: float = q_init[1]
                            c: float = q_init[2]
                            rho_0 *= Mixed2LevelState(i_qubit, a, b, c)
                        else:
                            raise Exception(
                                f"The initial state of site {i_qubit} is defined "
                                "using a tuple of an unsupported length."
                            )
                    else:
                        rho_0 *= get_operator_from_label(q_init, i_qubit)
                if h_x[i_qubit]:
                    H += (0.5 * h_x[i_qubit]) * Sx(i_qubit)
                if h_y[i_qubit]:
                    H += (0.5 * h_y[i_qubit]) * Sy(i_qubit)
                if h_z[i_qubit]:
                    H += (0.5 * h_z[i_qubit]) * Sz(i_qubit)
                if g_0[i_qubit]:
                    L_ops.append(Sp(i_qubit))
                    L_sig.append(Signal(g_0[i_qubit]))
                if g_1[i_qubit]:
                    L_ops.append(Sm(i_qubit))
                    L_sig.append(Signal(g_1[i_qubit]))
                if g_2[i_qubit]:
                    L_ops.append(Sz(i_qubit))
                    L_sig.append(Signal(g_2[i_qubit]))
                if g_3[i_qubit]:
                    L_ops.append(Sx(i_qubit))
                    L_sig.append(Signal(g_3[i_qubit]))
                if g_4[i_qubit]:
                    L_ops.append(Sy(i_qubit))
                    L_sig.append(Signal(g_4[i_qubit]))

            for i in r_qubits:
                for j in r_qubits:
                    if J[i, j]:
                        H += 0.5 * J[i, j] * (Sx(i) * Sx(j) + Sy(i) * Sy(j))
                    if J_z[i, j]:
                        H += 0.5 * J_z[i, j] * (Sz(i) * Sz(j))
            for i_qubit in _1q_indices:
                for s_op in _1q_components:
                    obs_1q.append(get_operator_from_label(s_op, i_qubit))
                    obs_1q_key.append((s_op, i_qubit))
            for _2q in _2q_indices:
                for s_op in _2q_components:
                    obs_2q.append(
                        get_operator_from_label(s_op[0], _2q[0])
                        * get_operator_from_label(s_op[1], _2q[1])
                    )
                    obs_2q_key.append((s_op, _2q[0], _2q[1]))
            for _3q in _3q_indices:
                for s_op in _3q_components:
                    obs_3q.append(
                        get_operator_from_label(s_op[0], _3q[0])
                        * get_operator_from_label(s_op[1], _3q[1])
                        * get_operator_from_label(s_op[2], _3q[2])
                    )
                    obs_3q_key.append((s_op, _3q[0], _3q[1], _3q[2]))
            for obs in custom_observables:
                gs_mat = build_matrices(gs, subsystem_dims)
                gate_matrix = build_matrices(id_op, subsystem_dims)
                if obs[0][1] == "g":
                    for obs_gate in obs[1]:  # Matrices are reversed w.r.t gates
                        gate_name = obs_gate[0]
                        if gate_name == "cz":
                            i = obs_gate[1]
                            j = obs_gate[2]
                            gate = 0.5 * (Sz(i) + Sz(j) - Sz(i) * Sz(j) + Id(i) * Id(j))
                        else:
                            gate = get_operator_from_label(gate_name, obs_gate[1])
                        gate_matrix = build_matrices(gate, subsystem_dims) @ gate_matrix
                    obs_mat = gate_matrix @ gs_mat @ gate_matrix.conj().T
                    obs_cu_mat.append(obs_mat)
                    obs_cu_key.append(obs[0][0])

            H_matrix = build_matrices(H, subsystem_dims)
            L_matrices = build_matrices(L_ops, subsystem_dims)
            hamiltonian = HamiltonianModel(operators=[H_matrix], signals=[Signal(1.0)])
            lindbladian = LindbladModel.from_hamiltonian(
                hamiltonian=hamiltonian,
                dissipator_operators=L_matrices,
                dissipator_signals=L_sig,
            )
            if s_load_files_prefix != "":  # Load initial state
                rho_0_mat = np.load(s_load_files_prefix + ".state.npy")
            else:
                rho_0_mat = build_matrices(rho_0, subsystem_dims)
                if init_cz_gates is not None:
                    # Apply CZ to all qubit tuples in the list
                    for q_pair in init_cz_gates:
                        i = q_pair[0]
                        j = q_pair[1]
                        CZ = 0.5 * (Sz(i) + Sz(j) - Sz(i) * Sz(j) + Id(i) * Id(j))
                        CZ_matrix = build_matrices(CZ, subsystem_dims)
                        rho_0_mat = CZ_matrix @ rho_0_mat @ CZ_matrix
                        # No dagger required as CZ is real diagonal.

            obs_1q_mat = build_matrices(obs_1q, subsystem_dims)
            obs_2q_mat = build_matrices(obs_2q, subsystem_dims)
            obs_3q_mat = build_matrices(obs_3q, subsystem_dims)
            y0 = rho_0_mat / np.trace(rho_0_mat)
            t_eval = np.arange(t_init, t_final, tau)
            if t_final not in t_eval:
                t_eval = np.concatenate((t_eval, [t_final]))
            sol = solve_lmde(
                lindbladian,
                t_span=Array((t_init, t_final)),
                y0=y0,
                t_eval=Array(t_eval),
                method=method,
                atol=atol,
                rtol=rtol,
            )

            if parameters.get("b_save_final_state", False):  # Write final state
                np.save(self.s_output_path + ".state", sol.y[-1])
            cut_off_observable = parameters.get("cut_off_observable", 0.0)

            file_1q = open(self.s_output_path + ".obs-1q.dat", "w")
            file_2q = open(self.s_output_path + ".obs-2q.dat", "w")
            file_3q = open(self.s_output_path + ".obs-3q.dat", "w")
            file_cu = open(self.s_output_path + ".obs-cu.dat", "w")
            file_gl = open(self.s_output_path + ".global.dat", "w")
            file_1q.write("#time\toperator\tindex\tvalue\n")
            file_2q.write("#time\toperator\tindex_1\tindex_2\tvalue\n")
            file_3q.write("#time\toperator\tindex_1\tindex_2\tindex_3\tvalue\n")
            file_cu.write("#time\tobservable\tvalue\n")
            file_gl.write("#time\tquantity\tvalue\n")

            n_times = len(sol.y)
            for t_i, sol_t in enumerate(sol.y):
                rho_t = DensityMatrix(sol_t)
                t = t_eval[t_i]

                for i_obs, _ in enumerate(obs_1q):
                    val = rho_t.expectation_value(obs_1q_mat[i_obs])
                    key = obs_1q_key[i_obs]
                    if abs(val.imag) > self.IMAGINARY_THRESHOLD:
                        self._print(
                            f"Warning: imaginary component {val.imag} in observable "
                            f"{key[0].upper()}, qubit {key[1] + 1}.\n"
                        )
                    if cut_off_observable and (abs(val) < cut_off_observable):
                        val = 0.0
                    file_1q.write(f"{t}\t{key[0].upper()}\t{key[1] + 1}\t{val.real}\n")
                file_1q.write("\n")
                file_1q.flush()

                for i_obs, _ in enumerate(obs_2q):
                    val = rho_t.expectation_value(obs_2q_mat[i_obs])
                    key = obs_2q_key[i_obs]
                    if abs(val.imag) > self.IMAGINARY_THRESHOLD:
                        self._print(
                            f"Warning: imaginary component {val.imag} in observable "
                            f"{key[0].upper()}, qubits ({key[1] + 1}, {key[2] + 1}).\n"
                        )
                    if cut_off_observable and (abs(val) < cut_off_observable):
                        val = 0.0
                    file_2q.write(
                        f"{t}\t{key[0].upper()}\t{key[1] + 1}\t{key[2] + 1}\t{val.real}\n"
                    )
                file_2q.write("\n")
                file_2q.flush()

                for i_obs, _ in enumerate(obs_3q):
                    val = rho_t.expectation_value(obs_3q_mat[i_obs])
                    key = obs_3q_key[i_obs]
                    if abs(val.imag) > self.IMAGINARY_THRESHOLD:
                        self._print(
                            f"Warning: imaginary component {val.imag} in observable "
                            f"{key[0].upper()}, qubits ({key[1] + 1}, {key[2] + 1}, "
                            f"{key[3] + 1}).\n"
                        )
                    if cut_off_observable and (abs(val) < cut_off_observable):
                        val = 0.0
                    file_3q.write(
                        f"{t}\t{key[0].upper()}\t{key[1] + 1}\t{key[2] + 1}\t{key[3] + 1}\t{val.real}\n"
                    )
                file_3q.write("\n")
                file_3q.flush()

                for obs_cu_matrix, obs_cu_name in zip(obs_cu_mat, obs_cu_key):
                    val = rho_t.expectation_value(obs_cu_matrix)
                    if abs(val.imag) > self.IMAGINARY_THRESHOLD:
                        self._print(
                            f"Warning: imaginary component {val.imag} in observable "
                            f"{obs_cu_name}.\n"
                        )
                    if cut_off_observable and (abs(val) < cut_off_observable):
                        val = 0.0
                    file_cu.write(f"{t}\t{obs_cu_name}\t{val.real}\n")
                    # print(f"{t}\t{obs_cu_name}\t{val.real}")
                file_cu.write("\n")
                file_cu.flush()

                rho = rho_t.data
                tr_rho = np.trace(rho).real
                if abs(tr_rho - 1) > self.TRACE_RHO_WARN_THRESHOLD:
                    self._print(f"Warning: Tr[rho] != 1 : {tr_rho}.\n")
                S_2 = -math.log(sum(scipy.linalg.eigvals(rho) ** 2).real)
                OSEE = np.nan
                bd = np.nan
                if t_i == n_times - 1:
                    ts1 = time.time()
                    duration = ts1 - ts0
                    self._print(f"Execution time: {duration}s")
                else:
                    duration = np.nan
                file_gl.write(f"{t}\ttr_rho\t{tr_rho}\n")
                file_gl.write(f"{t}\tS_2\t{S_2}\n")
                file_gl.write(f"{t}\tOSEE_center\t{OSEE}\n")
                file_gl.write(f"{t}\tmax_bond_dim\t{bd}\n")
                file_gl.write(f"{t}\tduration_ms\t{duration * 1000}\n")
                file_gl.write(f"\n")
                file_gl.flush()

            self.result = self.load_output(self.s_output_path)
        except Exception as e:
            self._print(str(e))
            raise e
        finally:
            self._close_file(self._log_file)
            self._log_file = None
            self._close_file(file_1q)
            self._close_file(file_2q)
            self._close_file(file_3q)
            self._close_file(file_cu)
            self._close_file(file_gl)

    def _get_parameter(
        self, s_key: str
    ) -> Union[None, Sized, Iterable, str, float, list, np.ndarray]:
        val = self.parameters.get(s_key, None)
        def_val = self.DEFAULT_PARAMETERS.get(s_key, None)

        if val is None and def_val is not None:
            val = def_val[0]
        if (
            (isinstance(val, str) or not isinstance(val, Container))
            and not isinstance(val, np.ndarray)
            and def_val is not None
            and def_val[1] != "s"
        ):
            # Not a scalar
            n_qubits = self.parameters["N"]
            if def_val[1] == "v":
                val = [val] * n_qubits  # create a list of identical values
            elif def_val[1] == "m":
                val = np.full((n_qubits, n_qubits), val)
                # create a matrix of identical values
        return val

    @staticmethod
    def _close_file(file):
        if file is not None:
            file.close()

    def _print(self, s_output: str):
        print(s_output)
        if self._log_file is not None:
            self._log_file.write(s_output + "\n")
            self._log_file.flush()

    def _virtual_verify_parameters(self, ignore_params: Optional[list] = None) -> str:
        """Returns a detailed Error message if parameters are not in the correct format.

        Args:
                ignore_params: A list with parameter names that this solver does not recognize, but
                        should be ignored in the verification (so that an error message for unknown parameters
                        is not issued). This is useful for derived classes.
        Returns:
                A detailed error message if parameters arguments are not in the correct format (which
                        is stated in the spec of the simulator). Otherwise, returns "" (checks passed).
        """
        return LindbladMatrixSolver.verify_parameters(self.parameters, ignore_params)

    @staticmethod
    def verify_parameters(
        parameters: dict, ignore_params: Optional[list] = None
    ) -> str:
        """Returns a detailed Error message if parameters are not in the correct format.

        Args:
                parameters: A dictionary of solver parameters.
                ignore_params: A list with parameter names that this solver does not recognize, but
                        should be ignored in the verification (so that an error message for unknown parameters
                        is not issued). Passed on to the base class (LindbladMPOSolver).
        Returns:
                A detailed error message if parameters are not in the correct format.
                Otherwise, returns "" (checks passed).
        """
        check_msg = ""
        unsupported_mpo_params = [
            "l_x",
            "l_y",
            "b_periodic_x",
            "b_periodic_y",
            "b_force_rho_trace",
            "force_rho_hermitian_step",
            "force_rho_hermitian_gates",
            "trotter_order",
            "b_initial_rho_compression",
            "b_apply_gate_compression",
            "b_collapse_equal_mixture",
            "max_dim_rho",
            "cut_off",
            "cut_off_rho",
        ]

        for key in dict.keys(parameters):
            if isinstance(parameters[key], str) and "" == parameters[key]:
                # ignore empty entrances/space holders <"">
                continue
            if key in unsupported_mpo_params:
                check_msg += (
                    "LindbladMatrixSolver Error 1000: "
                    + key
                    + " parameter is unsupported\n"
                )
                continue
            if key == "output_step":
                if not LindbladMPOSolver._is_int(parameters[key]):
                    check_msg += (
                        "LindbladMatrixSolver Error 1010: "
                        + key
                        + " should be an integer\n"
                    )
                    continue
                if parameters[key] < 0:
                    check_msg += (
                        "LindbladMatrixSolver Error 1020: "
                        + key
                        + " should be bigger/equal to 0 (integer)\n"
                    )
                    continue
            elif (key == "atol") or (key == "rtol"):
                if not LindbladMPOSolver.is_float(parameters[key]):
                    check_msg += (
                        "LindbladMatrixSolver Error 1030: " + key + " is not a float\n"
                    )
                    continue

        if ignore_params is None:
            ignore_params = []
        ignore_params.extend(["atol", "rtol", "method"])
        check_msg += LindbladMPOSolver.verify_parameters(parameters, ignore_params)
        return check_msg
