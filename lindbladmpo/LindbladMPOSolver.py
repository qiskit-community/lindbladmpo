# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Defines the main class of the packages, implementing the interface with the solver.
"""

import collections
import subprocess
import uuid
from math import isfinite
from typing import Dict, Optional
import platform
import os
import numpy as np


class LindbladMPOSolver:
    """Evolve multi-qubit Lindblad dynamics with a high-performance matrix-product-operators solver."""

    DEFAULT_EXECUTABLE_PATH = "/../bin/lindbladmpo"
    """Name and location of the solver executable. Will be appended with '.exe' for Windows."""

    DEFAULT_CYGWIN_PATH = "C:/cygwin64/bin/bash.exe"
    """Default path for the cygwin executable, which is used to invoke the solver (Windows only)."""

    def __init__(
        self,
        parameters: Optional[dict] = None,
        s_cygwin_path: Optional[str] = None,
        s_solver_path: Optional[str] = None,
    ):
        """Initialize the instance.

        Args:
                parameters: The model parameters.
                s_cygwin_path: On Windows only, indicates the cygwin executable path. A default
                        location will be assigned if this argument is not pass.
                s_solver_path: Indicates the solver executable path. A default location will be
                        assigned if this argument is not passed.
        """
        self.parameters = parameters
        self.s_input_file = ""
        self.s_output_path = ""
        s_cygwin_path, s_solver_path = self.process_default_paths(
            s_cygwin_path, s_solver_path
        )
        self.s_cygwin_path = s_cygwin_path
        self.s_solver_path = s_solver_path
        self.s_id_suffix = ""
        self.result = {}

    def solve(self):
        """Solves the simulation and loads the result dictionaries."""
        if self.s_input_file == "":
            self.build()
        exit_code = LindbladMPOSolver.execute(
            self.s_cygwin_path, self.s_solver_path, self.s_input_file
        )
        if exit_code != 0:
            raise Exception("There was an error executing the solver.")
        self.result = LindbladMPOSolver.load_output(self.s_output_path)

    @staticmethod
    def process_default_paths(
        s_cygwin_path: Optional[str] = None, s_solver_path: Optional[str] = None
    ) -> (str, str):
        """Returns the proper default values for the cygwin path and solver path according to
                the system platform, for each of those parameter that is None, otherwise the parameter
                is returned unchanged.
        Args:
                s_cygwin_path: the path for the cygwin executable (on Windows).
                s_solver_path: the path for the solver executable.
        Returns:
                (s_cygwin_path, s_solver_path): Default values for the cygwin path and solver path
                according to the system platform.
        """
        if s_cygwin_path is None or s_solver_path is None:
            s_solver_path1 = os.path.dirname(os.path.abspath(__file__))
            s_system = platform.system().lower()
            if s_system == "windows":
                # On Windows we execute the solver using the cygwin bash, using the default path:
                s_cygwin_path1 = LindbladMPOSolver.DEFAULT_CYGWIN_PATH

                # s_solver_path should be of the form "/cygdrive/c/ ... ", and we use below a path
                # relative to the current file's path in the package
                s_solver_path1 = s_solver_path1.replace(":", "")
                s_solver_path1 = s_solver_path1.replace("\\", "/")
                s_solver_path1 = "/cygdrive/" + s_solver_path1
                s_solver_path1 += LindbladMPOSolver.DEFAULT_EXECUTABLE_PATH + ".exe"
            else:
                s_cygwin_path1 = ""
                s_solver_path1 += LindbladMPOSolver.DEFAULT_EXECUTABLE_PATH
            if s_cygwin_path is None:
                s_cygwin_path = s_cygwin_path1
            if s_solver_path is None:
                s_solver_path = s_solver_path1
        return s_cygwin_path, s_solver_path

    def build(self, parameters: Optional[dict] = None):
        """Write the parameters dictionary to the `.input.txt` file the solver expects.

        Initializes also the member fields:
                self.s_input_file: File name of solver input file,
                self.s_output_prefix: Prefix for solver output path,
                self.s_id_suffix: The unique id suffix (possibly empty),
                All fields are initialized based on the user settings in the parameters dictionary
                (or default values if not assigned).

        Args:
                parameters: The model parameters.

        Raises:
                Exception: If the build encounters some validation errors, those are given in
                    the error message.
        """
        if parameters is not None:
            self.parameters = parameters
        parameters = self.parameters
        check_params = self._virtual_verify_parameters()
        # check if there is a problem with the input, if "" returned there is no problem
        if check_params != "":
            raise Exception(check_params)
        # check if the user defined an input file name
        s_output_path: str = parameters.get("output_files_prefix", "")
        if s_output_path == "" or s_output_path[-1] in ("/", ".", "\\"):
            s_output_path += "lindblad"
        b_uuid = parameters.get("b_unique_id", False)
        if b_uuid:
            s_uuid = uuid.uuid4().hex
            print("Generating a unique id for this simulation: " + s_uuid)
            parameters["unique_id"] = s_uuid
        else:
            s_uuid = parameters.get("unique_id", "")
        s_id_suffix = ""
        if s_uuid != "":
            s_id_suffix = "." + s_uuid
        s_output_path += s_id_suffix
        s_input_file = s_output_path + ".input.txt"

        b_bond_indices = False
        first_bond_indices = []
        second_bond_indices = []
        interactions = []
        if "J" in parameters.keys():
            if isinstance(parameters["J"], np.ndarray):
                interactions.append("J")
        if "J_z" in parameters.keys():
            if isinstance(parameters["J_z"], np.ndarray):
                interactions.append("J_z")
        if len(interactions) == 2:
            if parameters["J"].shape == parameters["J_z"].shape:
                b_bond_indices = True
                for i in range(parameters["J"].shape[0]):
                    for j in range(parameters["J"].shape[1]):
                        if parameters["J"][i, j] != 0 or parameters["J_z"][i, j] != 0:
                            first_bond_indices.append(i + 1)
                            second_bond_indices.append(j + 1)
            else:
                raise Exception("J and J_z are not of the same size.")
        elif len(interactions) == 1:
            b_bond_indices = True
            for i in range(parameters[interactions[0]].shape[0]):
                for j in range(parameters[interactions[0]].shape[1]):
                    if parameters[interactions[0]][i, j] != 0:
                        first_bond_indices.append(i + 1)
                        second_bond_indices.append(j + 1)

        print("Creating solver input file:")
        s_input_file = s_input_file.replace("\\", "/")
        print(s_input_file)
        file = open(s_input_file, "w")
        for key in parameters.keys():
            if key == "b_unique_id" or parameters[key] is None:
                pass
            elif key == "output_files_prefix":
                file.write(key + " = " + s_output_path + "\n")
            elif (
                (key in ("J", "J_z"))
                and isinstance(parameters[key], np.ndarray)
                and len(parameters[key]) > 1
            ):
                # check if to create bond indices arrays
                not_first_value = False
                file.write(key + " = ")
                for i in range(len(first_bond_indices)):
                    if not_first_value:
                        file.write(",")
                    file.write(
                        str(
                            parameters[key][
                                first_bond_indices[i] - 1, second_bond_indices[i] - 1
                            ]
                        )
                    )
                    not_first_value = True
                file.write("\n")
            elif (
                key == "init_pauli_state"
                or key == "init_product_state"
                or key == "apply_gates"
                or key == "1q_components"
                or key == "2q_components"
                or key == "3q_components"
            ):
                val = parameters[key]
                if isinstance(val, (int, float, tuple, str)):
                    val_list = [val]
                else:
                    val_list = val
                n_indices = len(val_list)
                file.write(key + " = ")
                for i_op, op in enumerate(val_list):
                    b_tuple = isinstance(op, tuple)
                    if isinstance(op, str):
                        file.write(str(op).strip("'"))
                    elif isinstance(op, (int, float)) or (b_tuple and len(op) == 1):
                        if key == "init_product_state":
                            file.write("p ")
                        s_val = op[0] if b_tuple else op
                        file.write(str(s_val))
                    elif b_tuple:
                        if key == "init_product_state":
                            if len(op) == 2:
                                file.write("q " + str(op[0]) + " " + str(op[1]))
                            else:
                                file.write(
                                    "r "
                                    + str(op[0])
                                    + " "
                                    + str(op[1])
                                    + " "
                                    + str(op[2])
                                )
                        elif key == "apply_gates":
                            file.write(str(op[0]) + " " + str(op[1]))
                            for j_op in range(
                                2, len(op)
                            ):  # qubit indices, 1-based in the file.
                                file.write(" " + str(op[j_op] + 1))
                    if i_op != n_indices - 1:
                        file.write(",")
                file.write("\n")
            elif key == "custom_observables" or key == "collapse":
                file.write(key + " = ")
                observables: list = parameters[key]
                n_observables = len(observables)
                for i_obs, (obs_def, obs_components) in enumerate(observables):
                    file.write(obs_def[0])
                    file.write(" ")
                    file.write(obs_def[1])
                    file.write(":")
                    n_components = len(obs_components)
                    for i_component, obs_component in enumerate(obs_components):
                        n_elements = len(obs_component)
                        for i_element, element in enumerate(obs_component):
                            if i_element == 0:
                                file.write(element)
                            else:
                                file.write(
                                    str(element + 1)
                                )  # Qubit index for 'g' observables
                            if i_element < n_elements - 1:
                                file.write(" ")
                        if i_component < n_components - 1:
                            file.write(",")
                    if i_obs < n_observables - 1:
                        file.write(";")
                file.write("\n")
            elif isinstance(parameters[key], np.ndarray):
                file.write(key + " = ")
                for i in range(parameters[key].shape[0]):
                    file.write(str(parameters[key][i]))
                    if i + 1 != parameters[key].shape[0]:
                        file.write(",")
                file.write("\n")
            elif key == "1q_indices":
                n_indices = len(parameters[key])
                file.write(key + " = ")
                for i_site, site in enumerate(parameters[key]):
                    file.write(str(site + 1))
                    # +1 because Python indices are 0-based, while iTensor's are 1-based
                    if i_site != n_indices - 1:
                        file.write(",")
                file.write("\n")
            elif (
                key == "2q_indices"
                or key == "init_graph_state"
                or key == "init_cz_gates"
            ):
                file.write(key + " = ")
                n_tuples = len(parameters[key])
                for i_2q_tuple, _2q_tuple in enumerate(parameters[key]):
                    file.write(str(_2q_tuple[0] + 1) + "," + str(_2q_tuple[1] + 1))
                    # +1 because Python indices are 0-based, while iTensor's are 1-based
                    if i_2q_tuple != n_tuples - 1:
                        file.write(",")
                file.write("\n")
            elif key == "3q_indices":
                file.write(key + " = ")
                n_tuples = len(parameters[key])
                for i_3q_tuple, _3q_tuple in enumerate(parameters[key]):
                    file.write(
                        str(_3q_tuple[0] + 1)
                        + ","
                        + str(_3q_tuple[1] + 1)
                        + ","
                        + str(_3q_tuple[2] + 1)
                    )
                    # +1 because Python indices are 0-based, while iTensor's are 1-based
                    if i_3q_tuple != n_tuples - 1:
                        file.write(",")
                file.write("\n")
            else:
                file.write(key + " = " + str(parameters[key]).strip("[]") + "\n")
        if b_bond_indices:
            file.write(
                "first_bond_indices = "
                + str(first_bond_indices).strip("[]").replace(" ", "")
                + "\n"
            )
            file.write(
                "second_bond_indices = "
                + str(second_bond_indices).strip("[]").replace(" ", "")
                + "\n"
            )
        file.close()
        self.s_input_file = s_input_file
        self.s_output_path = s_output_path
        self.s_id_suffix = s_id_suffix

    @staticmethod
    def execute(s_cygwin_path=None, s_solver_path=None, s_input_file="") -> int:
        """Execute the simulation solver.

        Args:
                s_cygwin_path : the path of the cygwin bash terminal execution
                s_solver_path : the path of the simulator executable
                s_input_file : input file path, if empty then the simulator is run using default values.

        Returns:
                exit code : the exit code of the solver.
        """
        s_cygwin_path, s_solver_path = LindbladMPOSolver.process_default_paths(
            s_cygwin_path, s_solver_path
        )
        if s_cygwin_path:
            call_string = s_cygwin_path + ' --login -c "'
        else:
            call_string = ""
        call_string += s_solver_path
        if s_input_file:
            call_string += " input_file '" + str(s_input_file) + "'"
            if s_cygwin_path:
                call_string += '"'
        print("Executing solver with command:")
        print("\t" + call_string + "\n")

        process = subprocess.Popen(call_string, shell=True)
        exit_code = process.wait()
        print(f"Solver process terminated with exit code {exit_code}.\n")
        return exit_code

    @staticmethod
    def load_output(s_output_path: str):
        """Read the three solver output files and returns a dictionary with the results.
        Args:
                s_output_path : prefix of the output files path. To this string the corresponding file
                        endings according to each output type will be appended.
        Returns:
                result : A dictionary with three dictionaries storing the different output types.
        """
        result = {}
        s_output_types = ["obs-1q", "obs-2q", "obs-3q", "obs-cu", "global"]
        for s_output_type in s_output_types:
            result[s_output_type] = LindbladMPOSolver._read_data_file(
                s_output_path, s_output_type
            )
        return result

    @staticmethod
    def _read_data_file(s_output_path: str, s_output_type: str) -> Dict:
        """Reads one of the solver output files and returns a dictionary with the data.
        Args:
                s_output_path : prefix of the output files path. To this string the corresponding file
                        endings according to each output type will be appended.
                s_output_type : A string defining the observable type, one of the 1-qubit, 2-qubits,
                        or global observables.
        Returns:
                result : A dictionary with the result.
        """
        full_filename = s_output_path + f".{s_output_type}.dat"
        result = collections.OrderedDict()
        if os.path.isfile(full_filename):
            print("Loading output data file: " + full_filename)
            file = open(full_filename, "r")
            file.readline()
            for line in file:
                words = line.strip().split()
                if not words:
                    continue
                LindbladMPOSolver._read_data_line(s_output_type, words, result)
            file.close()
        else:
            print("Skipping non-existing file: " + full_filename)
        return result

    @staticmethod
    def _read_data_line(s_output_type: str, words: list, result: Dict):
        t = float(words[0])
        op = words[1]
        val = float(words[-1])
        if s_output_type == "obs-1q":
            q_index1 = int(words[2]) - 1
            # data files are storing 1-based indices because of iTensor, while we use 0-based indices
            q_indices = (q_index1,)
        elif s_output_type == "obs-2q":
            q_index1 = int(words[2]) - 1
            # data files are storing 1-based indices because of iTensor, while we use 0-based indices
            q_index2 = int(words[3]) - 1
            q_indices = (q_index1, q_index2)
        elif s_output_type == "obs-3q":
            q_index1 = int(words[2]) - 1
            # data files are storing 1-based indices because of iTensor, while we use 0-based indices
            q_index2 = int(words[3]) - 1
            q_index3 = int(words[4]) - 1
            q_indices = (q_index1, q_index2, q_index3)
        elif s_output_type in ["obs-cu", "global"]:
            q_indices = ()
        else:
            raise Exception(f"Unknown output type {s_output_type}.")
        # The result dictionary is indexed by a tuple, first entry is a name, second entry is
        # a tuple of qubit indices - 0 indices for the global data, 1 for 1Q observables, 2 for 2Q.
        obs_data = result.get((op.lower(), q_indices), None)
        # obs_data is a tuple, first entry is a list of times, second entry holds the values.
        if obs_data is None:
            obs_data = (list(), list())
            result[(op.lower(), q_indices)] = obs_data
        # TODO: optimize the list appends
        obs_data[0].append(t)
        obs_data[1].append(val)

    @staticmethod
    # checks if the value is int (for cleaner code)
    def _is_int(value):
        return isinstance(value, int)

    @staticmethod
    def is_float(value):
        """checks if the value is a float (for cleaner code).
        Args:
                value : value to test.
        Returns:
                True if value is a float or int.
        """
        # in python terms the value <4> is not float type, in the simulator context float
        # can also be a python int:
        return isinstance(value, (float, int))

    @staticmethod
    # returns the number of qubits based on the given parameters, returns -1 if found an error
    def _get_number_of_qubits(parameters: Dict) -> int:
        if "N" in parameters:
            if LindbladMPOSolver._is_int(parameters["N"]):
                return parameters["N"]
        return -1

    def _virtual_verify_parameters(self, ignore_params: Optional[list] = None) -> str:
        """An overridable function that verifies the parameters by calling verify_parameters().

        Args:
                ignore_params: A list with parameter names that this solver does not recognize, but
                        should be ignored in the verification (so that an error message for unknown
                        parameters is not issued). This is useful for derived classes.
        Returns:
                A detailed error message if parameters arguments are not in the correct format (which
                        is stated in the spec of the simulator). Otherwise, returns "" (checks passed).
        """
        return LindbladMPOSolver.verify_parameters(self.parameters, ignore_params)

    @staticmethod
    def verify_parameters(
        parameters: dict, ignore_params: Optional[list] = None
    ) -> str:
        """Returns a detailed Error message if parameters are not in the correct format.

        Args:
                parameters: A dictionary of solver parameters.
                ignore_params: A list with parameter names that this solver does not recognize, but
                        should be ignored in the verification (so that an error message for unknown
                        parameters is not issued). This is mostly useful for derived subclasses.
        Returns:
                A detailed error message if parameters are not in the correct format.
                Otherwise, returns "" (checks passed).
        """
        check_msg = ""
        if parameters is None:
            check_msg += "Error 100: The `parameters` dictionary must be assigned\n"
            return check_msg
        if (
            ("N" not in parameters)
            or ("t_final" not in parameters)
            or ("tau" not in parameters)
        ):
            check_msg += (
                "Error 110: N, t_final and tau must be defined as they do not have default "
                "values\n"
            )
            return check_msg
        N = LindbladMPOSolver._get_number_of_qubits(parameters)
        for key in dict.keys(parameters):
            if parameters[key] is None or (
                isinstance(parameters[key], str) and "" == parameters[key]
            ):  # ignore empty entrances/space holders <"">
                continue
            flag_continue = False

            if key == "N":
                if not LindbladMPOSolver._is_int(parameters[key]):
                    check_msg += "Error 120: " + key + " should be an integer\n"
                    continue
                if parameters[key] <= 0:
                    check_msg += (
                        "Error 130: " + key + " should be bigger/equal to 1 (integer)\n"
                    )
                    continue
            elif key == "t_init" or key == "t_final" or key == "tau":
                if not LindbladMPOSolver.is_float(parameters[key]):
                    check_msg += "Error 140: " + key + " is not a float\n"
                    continue
                if key == "tau" and parameters[key] <= 0:
                    check_msg += "Error 150: " + key + " must be larger than 0\n"
                    continue
                if key == "t_init" and parameters[key] > parameters["t_final"]:
                    check_msg += (
                        "Error 151: " + key + " must be equal or smaller than t_final\n"
                    )
                    continue
            elif (key == "l_x") or (key == "l_y"):
                if not LindbladMPOSolver._is_int(parameters[key]):
                    check_msg += "Error 160: " + key + " should be an integer\n"
                    continue
                if parameters[key] < 0:
                    check_msg += (
                        "Error 170: "
                        + key
                        + " should be equal to or larger than 0 (integer)\n"
                    )
                    continue
            elif (
                key == "output_step"
                or key == "force_rho_hermitian_step"
                or key == "force_rho_hermitian_gates"
            ):
                if not LindbladMPOSolver._is_int(parameters[key]):
                    check_msg += "Error 180: " + key + " should be an integer\n"
                    continue
                if parameters[key] < 0:
                    check_msg += (
                        "Error 190: "
                        + key
                        + " should be equal to or larger than 0 (integer)\n"
                    )
                    continue
            elif (
                (key == "h_x")
                or (key == "h_y")
                or (key == "h_z")
                or (key == "g_0")
                or (key == "g_1")
                or (key == "g_2")
            ):
                if LindbladMPOSolver.is_float(parameters[key]):
                    continue
                if N == -1:
                    check_msg += (
                        "Error 200: " + key + " could not be validated because 'N' "
                        "(or alternatively l_x, l_y) are not "
                        "defined properly\n "
                    )
                    continue
                if isinstance(parameters[key], list):
                    if len(parameters[key]) != N:
                        check_msg += (
                            "Error 210: " + key + " is not a float / N-length list / "
                            "numpy array (of floats)\n"
                        )
                        continue
                    for element in parameters[key]:
                        if not LindbladMPOSolver.is_float(element):
                            check_msg += (
                                "Error 220: " + key + "is not a float / N-length list "
                                "/ numpy array (of floats)\n "
                            )
                            flag_continue = True
                            break
                    if flag_continue:
                        continue
                elif isinstance(parameters[key], np.ndarray):
                    if (str((parameters[key]).dtype).find("int") == -1) and (
                        str((parameters[key]).dtype).find("float") == -1
                    ):
                        check_msg += (
                            "Error 230: " + key + " is not a float / N-length list / "
                            "numpy array (of floats)\n"
                        )
                        continue
                    if parameters[key].size == 1:
                        continue
                    if (parameters[key].shape[0] != N) or (
                        parameters[key].shape[0] != parameters[key].size
                    ):
                        check_msg += (
                            "Error 240: " + key + " is not a float / N-length list / "
                            "numpy array (of floats)\n"
                        )
                        continue
                else:
                    check_msg += (
                        "Error 250: " + key + " is not a float / N-length list / numpy "
                        "array (of floats)\n"
                    )
                    continue
            elif (key == "J_z") or (key == "J"):
                if LindbladMPOSolver.is_float(parameters[key]):
                    continue
                if N == -1:
                    check_msg += (
                        "Error 260: " + key + " could not be validated because 'N' "
                        "(or alternatively l_x, l_y) are not "
                        "defined properly\n"
                    )
                    continue
                if isinstance(parameters[key], list):
                    if len(parameters[key]) != N:
                        check_msg += (
                            "Error 270: "
                            + key
                            + " should be a constant, or a square matrix"
                            " (nested lists/np.array) of N^2 floats\n "
                        )
                        continue
                    for lst in parameters[key]:
                        if not isinstance(lst, list):
                            check_msg += (
                                "Error 280: "
                                + key
                                + "should be a constant, or a square "
                                "matrix (nested lists/np.array) of "
                                "floats with a size N^2\n "
                            )
                            flag_continue = True
                            break
                        if len(lst) != N:
                            check_msg += (
                                "Error 290: "
                                + key
                                + "should be a constant, or a square matrix (nested "
                                "lists/np.array) with N^2 floats\n"
                            )
                            flag_continue = True
                            break
                        for val in lst:
                            if not LindbladMPOSolver.is_float(val):
                                check_msg += (
                                    "Error 300: "
                                    + key
                                    + "should be a constant, or a square matrix (nested "
                                    "lists/np.array) in the size of number_of_qubits^2 "
                                    "of floats\n"
                                )
                                flag_continue = True
                                break
                        if flag_continue:
                            break
                    if flag_continue:
                        continue
                elif isinstance(parameters[key], np.ndarray):
                    if (str((parameters[key]).dtype).find("int") == -1) and (
                        str((parameters[key]).dtype).find("float") == -1
                    ):
                        check_msg += (
                            "Error 310: "
                            + key
                            + "should be a constant, or a square matrix (nested "
                            "lists/np.array) in the size of number_of_qubits^2 of "
                            "floats\n"
                        )
                        continue
                    if parameters[key].size == 1:
                        continue
                    if parameters[key].shape[0] != N:
                        check_msg += (
                            "Error 320: "
                            + key
                            + "should be a constant, or a square matrix (nested "
                            "lists/np.array) in the size of number_of_qubits^2 of "
                            "floats\n"
                        )
                        continue
                    if parameters[key].shape[0] ** 2 != parameters[key].size:
                        check_msg += (
                            "Error 330: "
                            + key
                            + "should be a constant, or a square matrix (nested "
                            "lists/np.array) in the size of number_of_qubits^2 of "
                            "floats\n"
                        )
                        continue
                else:
                    check_msg += (
                        "Error 340: "
                        + key
                        + " should be a constant, or a square matrix (nested "
                        "list/np.array) in the size of number_of_qubits^2 of floats\n"
                    )
                    continue
            elif (
                key == "apply_gates" or key == "custom_observables" or key == "collapse"
            ):
                if (
                    not isinstance(parameters[key], tuple)
                    and not isinstance(parameters[key], list)
                    and not isinstance(parameters[key], np.ndarray)
                ):
                    check_msg += (
                        "Error 345: "
                        + key
                        + " must be a tuple or a list/ array of tuples\n"
                    )
                    continue
                custom_list = (
                    [parameters[key]]
                    if isinstance(parameters[key], tuple)
                    else parameters[key]
                )
                if key == "apply_gates":
                    for g_tuple in custom_list:
                        tuple_len = len(g_tuple)
                        if tuple_len < 3 or tuple_len > 4:
                            check_msg += (
                                "Error 346: every member of "
                                + key
                                + " must be of 3 or 4 elements\n"
                            )
                            continue
                        if (
                            not LindbladMPOSolver.is_float(g_tuple[0])
                            or not isinstance(g_tuple[1], str)
                            or not LindbladMPOSolver._is_int(g_tuple[2])
                            or (
                                tuple_len > 3
                                and not LindbladMPOSolver._is_int(g_tuple[3])
                            )
                        ):
                            check_msg += (
                                "Error 347: each member of "
                                + key
                                + " must be a tuple of the form"
                                " (time, gate name, qubit, [qubit])\n"
                            )
                            continue
                else:  # Hence key == "custom_observables" or key == "collapse"
                    b_is_collapse = key == "collapse"
                    for g_tuple in custom_list:
                        tuple_len = len(g_tuple)
                        if (
                            tuple_len != 2
                            or not isinstance(g_tuple[0], tuple)
                            or len(g_tuple[0]) != 2
                            or not isinstance(g_tuple[1], list)
                        ):
                            check_msg += (
                                "Error 341: every member of "
                                + key
                                + " must be a 2-tuple of a 2-tuple and a list\n"
                            )
                            continue
                        obs_type = g_tuple[0][1]
                        if not isinstance(g_tuple[0][0], str):
                            check_msg += (
                                "Error 342: each member of the first element of"
                                + key
                                + " must be a tuple of the form"
                                " (obs_name, obs_type)\n"
                            )
                            continue
                        if b_is_collapse:
                            if obs_type != "o":
                                check_msg += (
                                    "Error 342: each member of the first element of"
                                    + key
                                    + " must be a tuple of the form"
                                    " (obs_name, obs_type), with obs_type being 'o' to indicate"
                                    " a 1Q operator expansion\n"
                                )
                                continue
                        elif obs_type != "g" and obs_type != "o":
                            check_msg += (
                                "Error 342: each member of the first element of"
                                + key
                                + " must be a tuple of the form (obs_name, obs_type),"
                                " with obs_type being either 'g' or 'o' to indicate"
                                " a gate-based observable or a 1Q operator expansion\n"
                            )
                            continue
                        for o_tuple in g_tuple[1]:
                            # tuple_len = len(o_tuple)
                            if obs_type == "g" and (
                                not isinstance(o_tuple[0], str)
                                or not LindbladMPOSolver._is_int(o_tuple[1])
                            ):
                                check_msg += (
                                    "Error 343: each member of gate-based component of"
                                    + key
                                    + " must be a tuple of the form"
                                    " (gate_name, q0, q1, ...)\n"
                                )
                                continue
                            if obs_type == "o" and (
                                not isinstance(o_tuple[0], str)
                                or not LindbladMPOSolver._is_int(o_tuple[1])
                            ):
                                check_msg += (
                                    "Error 343: each member of an operator-based component of"
                                    + key
                                    + " must be a tuple of the form"
                                    " (operator name, qubit)\n"
                                )
                                continue

            elif (key == "init_pauli_state") or (key == "init_product_state"):
                if (
                    not isinstance(parameters[key], (str, float, tuple))
                    and not isinstance(parameters[key], list)
                    and not isinstance(parameters[key], np.ndarray)
                ):
                    check_msg += (
                        "Error 350: "
                        + key
                        + " must be a string, float, tuple or a list/ array of strings/floats/tuples\n"
                    )
                    continue
                init_list = (
                    [parameters[key]]
                    if isinstance(parameters[key], (str, float, tuple))
                    else parameters[key]
                )
                for q_init in init_list:
                    if isinstance(q_init, (float, int)) or (
                        isinstance(q_init, tuple) and len(q_init) == 1
                    ):
                        val = q_init[0] if isinstance(q_init, tuple) else q_init
                        if (
                            not LindbladMPOSolver.is_float(val)
                            or not isfinite(val)
                            or val < 0.0
                            or val > 1.0
                        ):
                            check_msg += (
                                "Error 361: a float or a length-1 tuple member of "
                                + key
                                + " represents a probability and must be between 0 and 1\n"
                            )
                        continue
                    if isinstance(q_init, tuple):
                        for val in q_init:
                            if not LindbladMPOSolver.is_float(val) or not isfinite(val):
                                check_msg += (
                                    "Error 362: the values in a tuple member of "
                                    + key
                                    + " must be valid numbers\n"
                                )
                        if len(q_init) == 2:
                            if q_init[0] < 0 or q_init[0] > np.pi:
                                check_msg += (
                                    "Error 363: the first value in a length-2 tuple of "
                                    + key
                                    + " represents a polar angle and must be in the range 0 to pi\n"
                                )
                        elif len(q_init) == 3:
                            if not (
                                0 <= q_init[0] <= 1
                                and -1 <= q_init[1] <= 1
                                and -1 <= q_init[2] <= 1
                            ):
                                check_msg += (
                                    "Error 364: a tuple member of "
                                    + key
                                    + " with three elements must contain valid entries of a"
                                    " density matrix\n"
                                )
                        else:
                            check_msg += (
                                "Error 365: a tuple member of "
                                + key
                                + " must be of 1, 2, or 3 elements\n"
                            )
                        continue
                    if not isinstance(q_init, str):
                        check_msg += (
                            "Error 360: each member of "
                            + key
                            + " must be a string, a float, or a tuple\n"
                        )
                        continue
                    allowed_init = ["+x", "-x", "+y", "-y", "+z", "-z", "id"]
                    if q_init.lower() not in allowed_init:
                        check_msg += (
                            "Error 370: "
                            + key
                            + " can only be one of: +x, -x, +y, -y, +z, -z, id\n"
                        )
                        continue
            elif (
                (key == "b_periodic_x")
                or (key == "b_periodic_y")
                or (key == "b_force_rho_trace")
                or (key == "b_unique_id")
                or (key == "b_quiet")
                or (key == "b_save_final_state")
                or (key == "b_initial_rho_compression")
                or (key == "b_apply_gate_compression")
            ):
                if not isinstance(parameters[key], bool):
                    check_msg += (
                        "Error 390: " + key + " should be a boolean True or False\n"
                    )
                    continue
            elif key == "trotter_order":
                if not LindbladMPOSolver._is_int(parameters[key]):
                    check_msg += "Error 400: " + key + " should be 2, 3 or 4\n"
                    continue
                if (
                    (parameters[key] != 2)
                    and (parameters[key] != 3)
                    and (parameters[key] != 4)
                ):
                    check_msg += "Error 401: " + key + " should be 2, 3 or 4\n"
                    continue
            elif key == "max_dim_rho":  # int
                if (
                    not LindbladMPOSolver._is_int(parameters[key])
                    or parameters[key] < 0
                ):
                    check_msg += (
                        "Error 410: " + key + " must be a non-negative integer\n"
                    )
                    continue
            elif (key == "cut_off") or (key == "cut_off_rho"):
                if not LindbladMPOSolver.is_float(parameters[key]):
                    check_msg += "Error 420: " + key + " is not a float\n"
                    continue
            elif key == "metadata":
                if not isinstance(parameters[key], str):
                    check_msg += "Error 422: " + key + " is not a string\n"
                    continue
                if "\n" in parameters[key]:
                    check_msg += (
                        "Error 423: "
                        "The metadata string cannot contain the new line "
                        "character code ('\\n'). Please reformat the string\n"
                    )
                    continue
            elif key == "load_files_prefix" or key == "output_files_prefix":
                if not isinstance(parameters[key], str):
                    check_msg += "Error 425: " + key + " is not a string\n"
                    continue
            elif key == "1q_components":
                x_c = 0
                y_c = 0
                z_c = 0
                if not isinstance(parameters[key], list):
                    check_msg += (
                        "Error 430: "
                        + key
                        + " should be a list of size 1,2,3 with x,y,z\n"
                    )
                    continue
                if len(parameters[key]) > 3:
                    check_msg += (
                        "Error 440: "
                        + key
                        + " should be a list of size 1,2,3 with x,y,z\n"
                    )
                    continue
                for val in parameters[key]:
                    if not isinstance(val, str):
                        check_msg += (
                            "Error 441: " + key + " only takes x,y,z (or a subset)\n"
                        )
                        flag_continue = True
                        break
                    val = str.lower(val)
                    if val == "x":
                        x_c += 1
                    elif val == "y":
                        y_c += 1
                    elif val == "z":
                        z_c += 1
                    else:
                        check_msg += (
                            "Error 450: " + key + " only takes x,y,z (or a subset)\n"
                        )
                        flag_continue = True
                        break
                if flag_continue:
                    continue
                if (x_c > 1) or (y_c > 1) or (z_c > 1):
                    check_msg += (
                        "Error 460: " + key + " only takes x,y,z (or a subset)\n"
                    )
                    continue
            elif key == "1q_indices":
                if parameters[key] != "":
                    if not isinstance(parameters[key], list):
                        check_msg += (
                            "Error 470: "
                            + key
                            + " should be an integer list (1,2,3,4..)\n"
                        )
                        continue
                    if N == -1:
                        check_msg += (
                            "Error 480: " + key + "could not be validated because 'N'"
                            " (or alternatively l_x,"
                            " l_y) are not defined properly\n "
                        )
                        continue
                    for element in parameters[key]:
                        if not LindbladMPOSolver._is_int(element):
                            check_msg += (
                                "Error 490: "
                                + key
                                + " should be an integer list (1,2,3,4..)\n"
                            )
                            flag_continue = True
                            break
                        if element >= N:
                            check_msg += (
                                "Error 500: "
                                + key
                                + " should be an integer list listing "
                                "qubits, therefore integers in the "
                                "range 0 to N-1\n"
                            )
                            flag_continue = True
                            break
                    if flag_continue:
                        continue
                    if len(parameters[key]) > N:
                        check_msg += (
                            "Error 510: "
                            + key
                            + " 's length should be equal/smaller than "
                            "the amount of qubits\n "
                        )
                        continue
                    if not len(set(parameters[key])) == len(parameters[key]):
                        check_msg += (
                            "Error 520: "
                            + key
                            + " 's List does not contain unique elements"
                        )
                        continue
            elif key == "2q_components":
                if not isinstance(parameters[key], list):
                    check_msg += (
                        "Error 530: "
                        + key
                        + "only accepts xx,yy,zz,xy,xz,yz (or a subset) "
                        "as a strings list\n"
                    )
                    continue
                if len(parameters[key]) > 6:
                    check_msg += (
                        "Error 540: "
                        + key
                        + " only accepts xx,yy,zz,xy,xz,yz (or a subset)\n"
                    )
                    continue
                check_me = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                for val in parameters[key]:
                    val = str.lower(val)
                    if val == "xx":
                        check_me[0] += 1
                    elif val == "yy":
                        check_me[1] += 1
                    elif val == "zz":
                        check_me[2] += 1
                    elif val == "xy":
                        check_me[3] += 1
                    elif val == "xz":
                        check_me[4] += 1
                    elif val == "yz":
                        check_me[5] += 1
                    elif val == "yx":
                        check_me[6] += 1
                    elif val == "zx":
                        check_me[7] += 1
                    elif val == "zy":
                        check_me[8] += 1
                    else:
                        check_msg += (
                            "Error 550: "
                            + key
                            + " only accepts string from xx, yy, zz, xy, "
                            "xz, yz (or a permutation thereof)\n"
                        )
                        flag_continue = True
                        break
                if flag_continue:
                    continue
                for check_val in check_me:
                    if check_val > 1:
                        check_msg += (
                            "Error 550: "
                            + key
                            + " only accepts strings from xx, yy, zz, xy, "
                            "xz, yz (or a permutation thereof)\n"
                        )
                        flag_continue = True
                        break
                if flag_continue:
                    continue
            elif key == "3q_components":
                if not isinstance(parameters[key], list):
                    check_msg += (
                        "Error 530: "
                        + key
                        + "only accepts xx,yy,zz,xy,xz,yz (or a subset) "
                        "as a strings list\n"
                    )
                    continue
                for val in parameters[key]:
                    val = str.lower(val)
                    b_ok = True
                    if len(val) != 3:
                        b_ok = False
                    else:
                        for c in val:
                            if c not in "xyz":
                                b_ok = False
                                break
                    if not b_ok:
                        check_msg += (
                            "Error 531: "
                            + key
                            + "only accepts length-3 Pauli strings\n"
                        )
                        flag_continue = True
                        break
                if flag_continue:
                    continue
            elif (
                key == "2q_indices"
                or key == "3q_indices"
                or key == "init_graph_state"
                or key == "init_cz_gates"
            ):  # expecting an integer tuples list
                if not isinstance(parameters[key], list):
                    check_msg += (
                        "Error 570: " + key + " should be a list of tuples"
                        " containing integers\n"
                    )
                    continue
                if N == -1:
                    check_msg += (
                        "Error 580: " + key + " could not be validated because 'N' "
                        "(or alternatively l_x, "
                        "l_y) are not defined properly\n"
                    )
                    continue
                tup_len = 3 if key == "3q_indices" else 2
                for tup in parameters[key]:
                    if not isinstance(tup, tuple):
                        check_msg += (
                            "Error 590: "
                            + key
                            + " should be a list of tuples containing integers\n "
                        )
                        flag_continue = True
                        break
                    if (
                        (len(tup) != tup_len)
                        or (not LindbladMPOSolver._is_int(tup[0]))
                        or (not LindbladMPOSolver._is_int(tup[1]))
                        or (tup_len == 3 and not LindbladMPOSolver._is_int(tup[2]))
                    ):
                        check_msg += (
                            "Error 600: "
                            + key
                            + f" should be a list of tuples of size {tup_len}, "
                            "containing integers\n "
                        )
                        flag_continue = True
                        break
                    if (tup[0] >= N) or (tup[1] >= N) or (tup_len == 3 and tup[2] >= N):
                        check_msg += (
                            "Error 610: "
                            + key
                            + f" should be a list of tuples of size {tup_len}, "
                            "containing integers equal to or smaller than "
                            "the total number of qubits\n "
                        )
                        flag_continue = True
                        break
                if flag_continue:
                    continue
                if not len(set(parameters[key])) == len(parameters[key]):
                    check_msg += "Error 630: " + key + " contains duplicate elements\n"
                    continue
            elif ignore_params is None or key not in ignore_params:
                check_msg += "Error: unknown parameter key passed: " + key + "\n"
        # End of: "for key in dict.keys(parameters)"

        return check_msg
