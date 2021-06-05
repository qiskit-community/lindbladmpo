import subprocess
import numpy as np


class MPOLindbladSolver:
    def __init__(self, parameters, s_cygwin_path, s_simulator_path):
        self.parameters = parameters
        if "output_file" in self.parameters:
            self.s_output_file = parameters["output_file"]
        else:
            self.s_output_file = "out"
        if "input_file" in self.parameters:
            self.s_input_file = parameters["input_file"]
        else:
            self.s_input_file = "input_file.txt"
        self.s_cygwin_path = s_cygwin_path
        self.s_simulator_path = s_simulator_path
        self.result = {}

    @staticmethod
    def _read_1q_output_into_dict(filename):
        """ Reads the 1 qubit output file and returns a dictionary with the values in the format:
        key (tuple) = (qubit (int), Axis (string), time (float)) : value = expectation value (float)
        Args:
            filename (string): file location
        Returns:
            result (dictionary): the arguments for the simulator
		"""
        full_filename = filename + ".1q_obs.dat"
        print("Loading 1-qubit observables from file:")
        print(full_filename)
        file = open(full_filename, "r")
        result = {}
        file.readline()
        for line in file:
            words = line.strip().split()
            if not words:
                continue
            words[0] = words[0].lstrip('#')
            result[(int(words[2]), words[0], float(words[1]))] = float(str(words[3]))
        file.close()
        return result

    @staticmethod
    def _read_2q_output_into_dict(filename):
        """ Reads the 2 qubit output file and returns a dictionary with the values in the format:
        key = (qubit1 (int), qubit2 (int), Axis1 (string), Axis2 (string) ,time (float)), value = expectation value (float)
        Args:
            filename (string): file location
        Returns:
            result (dictionary): the arguments for the simulator
		"""
        full_filename = filename + ".2q_obs.dat"
        print("Loading 2-qubit observables from file:")
        print(full_filename)
        file = open(full_filename, "r")
        result = {}
        file.readline()
        for line in file:
            words = line.strip().split()
            if not words:
                continue
            result[(int(words[2]), int(words[3]), words[0], float(words[1]))] = float(words[4])
        file.close()
        return result

    @staticmethod
    # checks if the value is int (for cleaner code)
    def _is_int(value):
        return isinstance(value, int)

    @staticmethod
    # checks if the value is a float (for cleaner code)
    def _is_float(value):
        # in python terms the value <4> is not float type, in the simulator context float can also be a python int:
        return isinstance(value, (float, int))

    @staticmethod
    # returns the number of qubits based on the given parameters, returns -1 if found an error
    def _get_number_of_qubits(dict_in):
        if ("l_x" in dict_in) and ("l_y" in dict_in):
            if MPOLindbladSolver._is_int(dict_in["l_x"]) and MPOLindbladSolver._is_int(dict_in["l_y"]):
                return dict_in["l_x"] * dict_in["l_y"]
        elif "N" in dict_in:
            if MPOLindbladSolver._is_int(dict_in["N"]):
                return dict_in["N"]
        return -1

    @staticmethod
    def _check_argument_correctness(dict_in):
        """Returns a detailed Error message if dict_in arguments are not in the correct format.
        Args:
            dict_in: a dictionary with the parameters for the simulation, keys as arguments and values as values.
        Returns:
            A detailed Error message if dict_in arguments are not in the correct format (which is stated in the spec of the simulator).
            else, returns "" (checks passed).
        Raises:
            Nothing, all issues will come out in the output returned.
		"""
        check_msg = ""
        # fixme for a better checking, assign the default values to any missing key \
        # as in the simulator it will be assigned anyway.. just remember to handle l_x l_y and N correctly (if N exists lx lx ignored)\
        # template: \
        # input_dict = {'N': 5, 't_final': 33, 'tau': 0.1, 'h_x': 1.111, 'h_y': [2.222, 4, 2e-2, 777777, -0.3], 'h_z': 3.333,
        #                   'J_z': [[2, 5, 6, 7, 9],
        #                           [4.55, -4.1, 12, -33, 10],
        #                           [4.55, -1.1, 17, 0, 10],
        #                           [4.55, -4.1, 61, -33, 10],
        #                           [4.55, -1.1, 11, -33, 10]], 'J': [[2, 5, 6, 7, 9],
        #                                                             [4.55, -4.1, 11, -33, 11e-9],
        #                                                             [4.55, -2.1, 11, 0, 10],
        #                                                             [6.25, -3.1, 2e-5, -33, 10],
        #                                                             [4.15, -4.1, 11, -33, 11]],
        #                   'g_0': [2.222, 4, 2e-2, 777777, -0.3], 'g_1': -99, 'g_2': -4e-10, 'init_Pauli_state': "+a", 'l_x': 5,
        #                   'l_y': 1, 'b_periodic_x': True, 'b_periodic_y': False, 'trotter_order': 3, 'max_dim': 100,
        #                   'max_dim_rho': 200, 'cut_off': 1e-10, 'cut_off_rho': 1e-11, 'b_force_rho_trace': False,
        #                   'b_force_rho_hermitian': True, 'output_step': 1, 'save_state_file': "C:\\output_of_my_simulation",
        #                   '1q_components': ["x", "y", "z"], '1q_indices': [1, 3, 4],
        #                   '2q_components': ["XX", "XY", "XZ", "YY", "YZ", "ZZ"], '2q_indices': [(1, 2), (3, 1), (2, 4), (3, 4)]}

        for key in dict.keys(dict_in):
            if "" == dict_in[key]:  # ignore empty entrances
                continue
            flag_continue = False

            if key == "N":
                if not MPOLindbladSolver._is_int(dict_in[key]):
                    check_msg += "Error: " + key + " should be an integer\n"
                    continue
                if dict_in[key] <= 0:
                    check_msg += "Error: " + key + " should be bigger/equal to 1 (integer)\n"
                    continue

            elif (key == "t_final") or key == "tau":
                if not MPOLindbladSolver._is_float(dict_in[key]):
                    check_msg += "Error: " + key + " is not a float\n"
                    continue
                if dict_in[key] <= 0:
                    check_msg += "Error: " + key + " must be bigger then 0\n"
                    continue

            elif (key == "l_x") or (key == "l_y"):
                if not MPOLindbladSolver._is_int(dict_in[key]):
                    check_msg += "Error: " + key + " should be an integer\n"
                    continue
                if dict_in[key] < 0:
                    check_msg += "Error: " + key + " should be bigger/equal to 1 (integer)\n"
                    continue

            elif key == "output_step":
                if not MPOLindbladSolver._is_int(dict_in[key]):
                    check_msg += "Error: " + key + " should be an integer\n"
                    continue
                if dict_in[key] < 0:
                    check_msg += "Error: " + key + " should be bigger/equal to 0 (integer)\n"
                    continue

            elif (key == "h_x") or (key == "h_y") or (key == "h_z") or \
                    (key == "g_0") or (key == "g_1") or (key == "g_2"):
                if MPOLindbladSolver._is_float(dict_in[key]):
                    continue
                number_of_qubits = MPOLindbladSolver._get_number_of_qubits(dict_in)
                if number_of_qubits == -1:
                    check_msg += "Error: " + key + " could not be validated because 'N' (or alternatively l_x, " \
                                                   "l_y) are not defined properly\n "
                    continue
                if isinstance(dict_in[key], list):
                    if len(dict_in[key]) != number_of_qubits:
                        check_msg += "Error: " + key + " is not a float / N size list / numpy array (of floats)\n"
                        continue
                    for element in dict_in[key]:
                        if not MPOLindbladSolver._is_float(element):
                            check_msg += "Error: " + key + " is not a float / N size list / numpy array (of floats)\n"
                            flag_continue = True
                            break
                    if flag_continue:
                        continue
                elif isinstance(dict_in[key], np.ndarray):
                    if (str((dict_in[key]).dtype).find("int") == -1) and (str((dict_in[key]).dtype).find("float") == -1):
                        check_msg += "Error: " + key + " is not a float / N size list / numpy array (of floats)\n"
                        continue
                    if dict_in[key].size == 1:
                        continue
                    if(dict_in[key].shape[0] != number_of_qubits) or (dict_in[key].shape[0] != dict_in[key].size):
                        check_msg += "Error: " + key + " is not a float / N size list / numpy array (of floats)\n"
                        continue
                else:
                    check_msg += "Error: " + key + " is not a float / N size list / numpy array (of floats)\n"
                    continue

            elif (key == "J_z") or (key == "J"):
                if MPOLindbladSolver._is_float(dict_in[key]):
                    continue
                number_of_qubits = MPOLindbladSolver._get_number_of_qubits(dict_in)
                if number_of_qubits == -1:
                    check_msg += "Error: " + key + " could not be validated because 'N' (or alternatively l_x, " \
                                                   "l_y) are not defined properly\n"
                    continue
                if isinstance(dict_in[key], list):
                    if len(dict_in[key]) != number_of_qubits:
                        check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                       "list/np.array) in the size of number_of_qubits^2 of floats\n"
                        continue
                    for lst in dict_in[key]:
                        if not isinstance(lst, list):
                            check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                        "list/np.array) in the size of number_of_qubits^2 of floats\n"
                            flag_continue = True
                            break
                        if len(lst) != number_of_qubits:
                            check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                        "list/np.array) in the size of number_of_qubits^2 of floats\n"
                            flag_continue = True
                            break
                        for val in lst:
                            if not MPOLindbladSolver._is_float(val):
                                check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                         "list/np.array) in the size of number_of_qubits^2 of floats\n"
                                flag_continue = True
                                break
                        if flag_continue:
                            break
                    if flag_continue:
                        continue
                elif isinstance(dict_in[key], np.ndarray):
                    if (str((dict_in[key]).dtype).find("int") == -1) and (str((dict_in[key]).dtype).find("float") == -1):
                        check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                       "list/np.array) in the size of number_of_qubits^2 of floats\n"
                        continue
                    if dict_in[key].size == 1:
                        continue
                    if dict_in[key].shape[0] != number_of_qubits:
                        check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                       "list/np.array) in the size of number_of_qubits^2 of floats\n"
                        continue
                    if dict_in[key].shape[0]**2 != dict_in[key].size:
                        check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                       "list/np.array) in the size of number_of_qubits^2 of floats\n"
                        continue
                else:
                    check_msg += "Error: " + key + " should be a constant, or a square matrix (nested " \
                                                   "list/np.array) in the size of number_of_qubits^2 of floats\n"
                    continue

            elif key == "init_Pauli_state":
                if not isinstance(dict_in[key], str):
                    check_msg += "Error: " + key + " is not a string\n"
                    continue
                str_as_lst = list(dict_in[key])
                if len(str_as_lst) != 2:
                    check_msg += "Error: " + key + " is a string but not of length 2 !\n"
                    continue
                if (str_as_lst[0] != '+') and (str_as_lst[0] != "-"):
                    check_msg += "Error: " + key + " can only be one of these: +x, -x, +y, -y, +z, -z\n"
                    continue
                if (str_as_lst[1] != "x") and (str_as_lst[1] != "y") and (str_as_lst[1] != "z"):
                    check_msg += "Error: " + key + " can only be one of these: +x, -x, +y, -y, +z, -z\n"
                    continue

            elif ((key == "b_periodic_x") or (key == "b_periodic_y") or (key == "b_force_rho_trace") or (
                    key == "b_force_rho_hermitian")):
                if not isinstance(dict_in[key], bool):
                    check_msg += "Error: " + key + " should be a boolean True or False as its a switch\n"
                    continue

            elif key == "trotter_order":
                if ((not MPOLindbladSolver._is_int(dict_in[key])) and dict_in[key] != 2 and dict_in[key] != 3 and
                        dict_in[key] != 4):
                    check_msg += "Error: " + key + " should be 2,3 or 4\n"
                    continue

            elif (key == "max_dim") or (key == "max_dim_rho"):  # int
                if not MPOLindbladSolver._is_int(dict_in[key]):
                    check_msg += "Error: " + key + " should be an integer\n"
                    continue

            elif (key == "cut_off") or (key == "cut_off_rho"):
                if not MPOLindbladSolver._is_float(dict_in[key]):
                    check_msg += "Error: " + key + " is not a small float format (ae-b where a and b are numbers)\n"
                    continue

            elif key == "save_state_file":
                pass  # fixme: add a check for windows/linux/mac path string validation
                      # allow only letters, numbers and underscore ? _   (should be enough)
            elif key == "output_file":
                pass
            elif key == "input_file":
                pass
            elif key == "1q_components":
                x_c = 0
                y_c = 0
                z_c = 0
                if not isinstance(dict_in[key], list):
                    check_msg += "Error: " + key + " should get a list of sizes 1,2,3 with x,y,z)\n"
                    continue
                if len(dict_in[key]) > 3:
                    check_msg += "Error: " + key + " should get a list of sizes 1,2,3 with x,y,z)\n"
                    continue
                for val in dict_in[key]:
                    val = str.lower(val)
                    if val == "x":
                        x_c += 1
                    elif val == "y":
                        y_c += 1
                    elif val == "z":
                        z_c += 1
                    else:
                        check_msg += "Error: " + key + " only gets x,y,z (or a subset)\n"
                        flag_continue = True
                        break
                if flag_continue:
                    continue
                if (x_c > 1) or (y_c > 1) or (z_c > 1):
                    check_msg += "Error: " + key + " only gets x,y,z (or a subset)\n"
                    continue

            elif key == "1q_indices":
                if dict_in[key] != "":
                    if not isinstance(dict_in[key], list):
                        check_msg += "Error: " + key + " should be an integer list (1,2,3,4..)\n"
                        continue
                    number_of_qubits = MPOLindbladSolver._get_number_of_qubits(dict_in)
                    if number_of_qubits == -1:
                        check_msg += "Error: " + key + "could not be validated because 'N' (or alternatively l_x, " \
                                                       "l_y) are not defined properly\n"
                        continue
                    for element in dict_in[key]:
                        if not MPOLindbladSolver._is_int(element):
                            check_msg += "Error: " + key + " should be an integer list (1,2,3,4..)\n"
                            flag_continue = True
                            break
                        if element >= number_of_qubits:
                            check_msg += "Error: " + key + " should be an integer list listing qubits, therefore " \
                                                           "number range is 0 to 'num of qubits'-1\n"
                            flag_continue = True
                            break
                    if flag_continue:
                        continue
                    if len(dict_in[key]) > number_of_qubits:
                        check_msg += "Error: " + key + " 's length should be smaller/equal then the amount of qubits\n"
                        continue
                    if not len(set(dict_in[key])) == len(dict_in[key]):
                        check_msg += "Error: " + key + " 's List does not contains all unique elements"
                        continue

            elif key == "2q_components":
                if not isinstance(dict_in[key], list):
                    check_msg += "Error: " + key + "only receives xx,yy,zz,xy,xz,yz (or a subset) as a strings " \
                                                   "list\n"
                    continue
                if len(dict_in[key]) > 6:
                    check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
                    continue
                check_me = [0, 0, 0, 0, 0, 0]
                for val in dict_in[key]:
                    val = str.lower(val)
                    if val == "xx":
                        check_me[0] += 1
                    elif val == "yy":
                        check_me[1] += 1
                    elif val == "zz":
                        check_me[2] += 1
                    elif (val == "xy") or (val == "yx"):
                        check_me[3] += 1
                    elif (val == "xz") or (val == "zx"):
                        check_me[4] += 1
                    elif (val == "yz") or (val == "zy"):
                        check_me[5] += 1
                    else:
                        check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
                        flag_continue = True
                        break
                if flag_continue:
                    continue
                for check_val in check_me:
                    if check_val > 1:
                        check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
                        flag_continue = True
                        break
                if flag_continue:
                    continue

            elif key == "2q_indices":  # expecting an integer tuples list
                if not isinstance(dict_in[key], list):
                    check_msg += "Error: " + key + " should be an list of tuples of size 2, containing integer\n"
                    continue
                number_of_qubits = MPOLindbladSolver._get_number_of_qubits(dict_in)
                if number_of_qubits == -1:
                    check_msg += "Error: " + key + " could not be validated because 'N' (or alternatively l_x, " \
                                                   "l_y) are not defined properly\n"
                    continue
                for tup in dict_in[key]:
                    if not isinstance(tup, tuple):
                        check_msg += "Error: " + key + " should be an list of tuples of size 2, containing " \
                                                       "integer\n "
                        flag_continue = True
                        break
                    if ((not MPOLindbladSolver._is_int(tup[0])) or (not MPOLindbladSolver._is_int(tup[1])) or (
                            len(tup) != 2)):
                        check_msg += "Error: " + key + " should be an list of tuples of size 2, containing " \
                                                       "integers\n "
                        flag_continue = True
                        break
                    if (tup[0] >= number_of_qubits) or (tup[1] >= number_of_qubits):
                        check_msg += "Error: " + key + " should be an list of tuples of size 2, containing " \
                                                       "integers smaller/equal then the total number of qubits\n "
                        flag_continue = True
                        break
                if flag_continue:
                    continue

                if len(dict_in[key]) > number_of_qubits ** 2:
                    check_msg += "Error: " + key + " 's length should be smaller then the amount of qubits^2\n"
                    continue
                if not len(set(dict_in[key])) == len(dict_in[key]):
                    check_msg += "Error: " + key + " 's List does not contains all unique elements"
                    continue

            else:
                check_msg += "Error: invalid key (parameter) inserted: " + key + "\n"
        # End of: "for key in dict.keys(dict_in)"
        # More cross-parameter checks:
        if ("t_final" in dict_in) and ("tau" in dict_in):
            if (MPOLindbladSolver._is_float(dict_in["tau"])) and (MPOLindbladSolver._is_float(dict_in["t_final"])):
                if (dict_in["tau"] > 0) and (dict_in["t_final"] > 0):
                    if dict_in["tau"] > dict_in["t_final"]:
                        check_msg += "Error: t_final (total time) is smaller then tau (time step for time evolution)\n"
                    elif "output_step" in dict_in:  # fixme: is this check relevant ? ( output_step*tau < t_final) if not, delete the whole elif
                        if MPOLindbladSolver._is_int(dict_in["output_step"]):
                            if dict_in["output_step"] > 0:
                                if dict_in["output_step"] * dict_in["tau"] > dict_in["t_final"]:
                                    check_msg += "Error: Output_step multiplied by tau, is bigger then t_final (output " \
                                                 "step in units of tau, times tau is bigger then the simulation " \
                                                 "time)\n "
        return check_msg

    @staticmethod
    def build_input_file(parameters):
        """ Writing the input parameters from the input dictionary to txt file
        Args:
            parameters (dictionary): the arguments for the simulator
        Returns:
		"""
        check_output = MPOLindbladSolver._check_argument_correctness(parameters)
        if check_output != "":
            print(check_output)
            raise Exception(
                check_output)  # fixme: is there a problem with this exception ? if raised, it stops the run for me and not generate an input file as intended (you said it wont stop the input file generation)
        if "input_file" in parameters.keys():
            file_name = parameters["input_file"]
        else:
            file_name = "input_file.txt"
        print("Preparing solver input file:")
        print(file_name)
        file = open(file_name, "w")
        for key in parameters.keys():
            file.write(key + " = " + str(parameters[key]).strip("[]") + "\n")
        file.close()

    @staticmethod
    def execute(s_cygwin_path, s_simulator_path, s_input_file=""):
        """ Execute the C++ simulator
        Args:
            s_cygwin_path (string): the address of the cygwin bash terminal execution
            s_simulator_path (string): the address of the simulator .exe
            s_input_file (string): input file location and name, if empty then the simulator is running on default values

        Returns:
		"""
        if s_cygwin_path:
            # call_string = 'cmd /k ' + s_cygwin_path + " --login -i -c \""
            call_string = s_cygwin_path + " --login -c \""
        else:
            call_string = ''
        call_string += s_simulator_path
        if s_input_file:
            call_string += " input_file " + str(s_input_file)
        print("Executing solver with command:")
        print("\t" + call_string + "\n")

        process = subprocess.Popen(call_string)
        exit_code = process.wait()
        print(f"Solver process terminated with exit code {exit_code}")

    @staticmethod
    def load_output(s_output_file):
        result = {'1q': MPOLindbladSolver._read_1q_output_into_dict(s_output_file),
                  '2q': MPOLindbladSolver._read_2q_output_into_dict(s_output_file)}
        return result

    def solve(self):
        self.build_input_file(self.parameters)
        self.execute(self.s_cygwin_path, self.s_simulator_path, self.s_input_file)
        self.result = self.load_output(self.s_output_file)
