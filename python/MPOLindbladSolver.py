import subprocess
import os
import signal
import time


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
	def read_1q_output_into_dict(filename):
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
	def read_2q_output_into_dict(filename):
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
	def is_int(value):
		return isinstance(value, int)

	@staticmethod
	def is_float(value):
		return (isinstance(value, float) or isinstance(value, int))

	@staticmethod
	def is_a_float_list(lst):
		if (not isinstance(lst, list)):
			return False
		for element in lst:
			if (not MPOLindbladSolver.is_float(element)):
				return False
		return True

	@staticmethod
	def get_number_of_qubits(dict_in):
		if (("l_x" in dict_in) and ("l_y" in dict_in)):
			if (MPOLindbladSolver.is_int(dict_in["l_x"]) and MPOLindbladSolver.is_int(dict_in["l_y"])):
				return (dict_in["l_x"] * dict_in["l_y"])
		elif ("n_qubits" in dict_in):
			if (MPOLindbladSolver.is_int(dict_in["n_qubits"])):
				return dict_in["n_qubits"]
		return -1

	@staticmethod
	def check_argument_correctness(dict_in):
		"""
		Returns a detailed Error message if dict_in agruments are not in the correct format.
			Args:
				dict_in: a dictionary with the parameters for the simulaton, keys as agruments and values as values.
			Returns:
				A detailed Error message if dict_in agruments are not in the correct format (which is stated in the spec of the simulator).
				else, returns "" (checks passed).
			Raises:
				Nothing, all issues will come out in the output returned.
		"""

		check_msg = ""
		for key in dict.keys(dict_in):

			if (key == "n_qubits"):
				if (not MPOLindbladSolver.is_int(dict_in[key])):
					check_msg += "Error: " + key + " should be an integer\n"
				else:
					if (dict_in[key] < 1):
						check_msg += "Error: " + key + " should be bigger/equal to 1 (integer)\n"

			elif ((key == "t_final") or (key == "output_step")):
				if (not MPOLindbladSolver.is_int(dict_in[key])):
					check_msg += "Error: " + key + " should be an integer\n"

			elif (key == "tau"):
				if (MPOLindbladSolver.is_float(dict_in[key])):
					if (dict_in[key] <= 0):
						check_msg += "Error: " + key + " must be bigger then 0\n"
				else:
					check_msg += "Error: " + key + " is not a float\n"

			elif ((key == "h_x") or (key == "h_y") or (key == "h_z") or (key == "g_0") or (key == "g_1") or (
					key == "g_2")):
				if (not MPOLindbladSolver.is_float(dict_in[key])):
					if (isinstance(dict_in[key], list)):
						number_of_qubits = MPOLindbladSolver.get_number_of_qubits(dict_in)
						if (number_of_qubits == -1):
							check_msg += "Error: " + key + " could not be validated because 'n_qubits' (or alternativly l_x, l_y) are not defined properly\n"
						else:
							if (len(dict_in[key]) != number_of_qubits):
								check_msg += "Error: " + key + " is not a float or a N_Qubits size list (of floats)\n"
							else:
								for element in dict_in[key]:
									if (not MPOLindbladSolver.is_float(element)):
										check_msg += "Error: " + key + " is not a float or a N_Qubits size list (of floats)\n"
										break
					else:
						check_msg += "Error: " + key + " is not a float or a N_Qubits size list (of floats)\n"

			elif ((key == "J_z") or (key == "J")):
				if ((dict_in[key] != "") and (not MPOLindbladSolver.is_float(dict_in[key]))):
					if (not isinstance(dict_in[key], list)):
						check_msg += "Error: " + key + " should be either empty/const/square matrix (floats)\n"
					else:
						number_of_qubits = MPOLindbladSolver.get_number_of_qubits(dict_in)
						if (number_of_qubits == -1):
							check_msg += "Error: " + key + " could not be validated because 'n_qubits' (or alternativly l_x, l_y) are not defined properly\n"
						else:
							if (len(dict_in[key]) != number_of_qubits):
								check_msg += "Error: " + key + " should be either empty/const/square matrix in the size of qubits^2 (floats)\n"
							else:
								for lst in dict_in[key]:
									if (not isinstance(lst, list)):
										check_msg += "Error: " + key + " should be either empty/const/square matrix (floats)\n"
										break
									elif (len(lst) != number_of_qubits):
										check_msg += "Error: " + key + " should be either empty/const/square matrix in the size of qubits^2 (floats)\n"
										break
									else:
										for val in lst:
											if (not MPOLindbladSolver.is_float(val)):
												check_msg += "Error: " + key + " should be either empty/const/square matrix (floats)\n"
												break

			elif (key == "init_pure_state"):
				if (not isinstance(dict_in[key], str)):
					check_msg += "Error: " + key + " is not a string\n"
				else:
					if (len(dict_in[key]) != 2):
						check_msg += "Error: " + key + " is a string but not of length 2 !\n"

			elif ((key == "b_periodic_x") or (key == "b_periodic_y") or (key == "b_force_rho_trace") or (
					key == "b_force_rho_hermitian")):
				if (not isinstance(dict_in[key], bool)):
					check_msg += "Error: " + key + " should be a boolian True or False as its a switch\n"

			elif (key == "trotter_order"):
				if ((not MPOLindbladSolver.is_int(dict_in[key])) and dict_in[key] != 2 and dict_in[key] != 3 and
						dict_in[key] != 4):
					check_msg += "Error: " + key + " should be 2,3 or 4\n"

			elif ((key == "max_dim") or (key == "max_dim_rho")):  # int
				if (not MPOLindbladSolver.is_int(dict_in[key])):
					check_msg += "Error: " + key + " should be an integer\n"

			elif ((key == "cut_off") or (key == "cut_off_rho")):
				if (not MPOLindbladSolver.is_float(dict_in[key])):
					check_msg += "Error: " + key + " is not a small float format (ae-b where a and b are numbers)\n"

			elif (key == "save_state_file"):
				pass  # fixme: add a check for windows/linux/mac path string validation
			elif (key == "output_file"):
				pass
			elif (key == "input_file"):
				pass
			elif (key == "1q_components"):
				if (dict_in[key] != ""):
					x_c = 0
					y_c = 0
					z_c = 0
					if (not isinstance(dict_in[key], list)):
						check_msg += "Error: " + key + " should get a of sizes 1,2,3 with x,y,z)\n"
					else:
						if (len(dict_in[key]) > 3):
							check_msg += "Error: " + key + " should get a of sizes 1,2,3 with x,y,z)\n"
						else:
							for val in dict_in[key]:
								val = str.lower(val)
								if ((val == "x") and (x_c != 1)):
									x_c += 1

								elif ((val == "y") and (y_c != 1)):
									y_c += 1

								elif ((val == "z") and (z_c != 1)):
									z_c += 1

								else:
									check_msg += "Error: " + key + " only gets x,y,z (or a subset)\n"
									break

			elif (key == "1q_sites"):
				if (dict_in[key] != ""):
					if (not isinstance(dict_in[key], list)):
						check_msg += "Error: " + key + " should be an integer list (1,2,3,4..)\n"
					else:
						for element in dict_in[key]:
							if (not MPOLindbladSolver.is_int(element)):
								check_msg += "Error: " + key + " should be an integer list (1,2,3,4..)\n"
								break

			elif (key == "2q_components"):
				if (dict_in[key] != ""):
					chkm = [0, 0, 0, 0, 0, 0]
					if (not isinstance(dict_in[key], list)):
						check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset) as a strings list\n"
					else:
						if (len(dict_in[key]) > 6):
							check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
						else:
							for val in dict_in[key]:
								val = str.lower(val)
								if ((val == "xx") and (chkm[0] != 1)):
									chkm[0] += 1
								elif ((val == "yy") and (chkm[1] != 1)):
									chkm[1] += 1
								elif ((val == "zz") and (chkm[2] != 1)):
									chkm[2] += 1
								elif (((val == "xy") or (val == "yx")) and (chkm[3] != 1)):
									chkm[3] += 1
								elif (((val == "xz") or (val == "zx")) and (chkm[4] != 1)):
									chkm[4] += 1
								elif (((val == "yz") or (val == "zy")) and (chkm[5] != 1)):
									chkm[5] += 1
								else:
									check_msg += "Error: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
									break

			elif (key == "2q_sites"):  # expecting an integer tuples list
				if (dict_in[key] != ""):
					if (not isinstance(dict_in[key], list)):
						check_msg += "Error: " + key + " should be an list of tuples of size 2, containing ingeter\n"
					else:
						for tup in dict_in[key]:
							if (not isinstance(tup, tuple)):
								check_msg += "Error: " + key + " should be an list of tuples of size 2, containing ingeter\n"
								break
							if ((not MPOLindbladSolver.is_int(tup[0])) or (not MPOLindbladSolver.is_int(tup[1])) or (
									len(tup) != 2)):
								check_msg += "Error: " + key + " should be an list of tuples of size 2, containing ingeters\n"
								break

			else:
				check_msg += "Error: " + "Error: user inserted invalid key (agrument): " + key + "\n"
		return check_msg

	@staticmethod
	def build_input_file(parameters):
		""" Writing the input parameters from the input dictionary to txt file
			Args:
				parameters (dictionary): the arguments for the simulator
			Returns:
		"""
		check_output = MPOLindbladSolver.check_argument_correctness(parameters)
		if check_output != "":
			print(check_output)
			return check_output
		if "input_file" in parameters.keys():
			file_name = parameters["input_file"]
		else:
			file_name = "input_file.txt"
		print("Preparing solver input file:")
		print(file_name)
		file = open(file_name, "w")
		for key in parameters.keys():
			file.write(key + " = " + str(parameters[key]) + "\n")
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
		result = {}
		result['1q'] = MPOLindbladSolver.read_1q_output_into_dict(s_output_file)
		result['2q'] = MPOLindbladSolver.read_2q_output_into_dict(s_output_file)
		return result

	def solve(self):
		self.build_input_file(self.parameters)
		self.execute(self.s_cygwin_path, self.s_simulator_path, self.s_input_file)
		self.result = self.load_output(self.s_output_file)
