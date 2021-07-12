# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import subprocess
import uuid
from typing import Dict
import numpy as np
import platform
import os


class LindbladMPOSolver:
	def __init__(self, parameters, s_cygwin_path = None, s_solver_path = None):
		self.parameters = parameters
		self.s_input_file = ''
		self.s_output_path = ''
		s_cygwin_path, s_solver_path = self.process_default_paths(s_cygwin_path, s_solver_path)
		self.s_cygwin_path = s_cygwin_path
		self.s_solver_path = s_solver_path
		self.s_id_suffix = ''
		self.s_N_suffix = ''
		self.result = {}

	def solve(self):
		try:
			self.s_input_file, self.s_output_path, self.s_id_suffix, self.s_N_suffix = self.build_input_file(self.parameters)
		except:
			raise Exception("There was an error creating the input file, aborting program")
		exit_code = self.execute(self.s_cygwin_path, self.s_solver_path, self.s_input_file)
		if exit_code != 0:
			raise Exception("There was an error executing the solver, aborting program.")
		self.result = self.load_output(self.s_output_path)

	@staticmethod
	def process_default_paths(s_cygwin_path = None, s_solver_path = None) -> (str, str):
		""" Returns the proper default values for the cygwin path and solver path according to the system
			platform, for each of those parameter that is None, otherwise the parameter is returned unchanged.
		Args:
			s_cygwin_path: the path for the cygwin executable (on Windows).
			s_solver_path: the path for the solver executable.
		Returns:
			(s_cygwin_path, s_solver_path): Default values for the cygwin path and solver path according to
			the system platform.
		"""
		if s_cygwin_path is None or s_solver_path is None:
			s_solver_path1 = os.path.dirname(os.path.abspath(__file__))
			s_system = platform.system().lower()
			if s_system == 'windows':
				# On Windows we execute the solver using the cygwin bash. The default path is the following.
				s_cygwin_path1 = "C:/cygwin64/bin/bash.exe"

				# s_solver_path should be of the form "/cygdrive/c/ ... ", and we use below a path relative
				# to the current file's path in the package
				s_solver_path1 = s_solver_path1.replace(':', '')
				s_solver_path1 = s_solver_path1.replace('\\', '/')
				s_solver_path1 = "/cygdrive/" + s_solver_path1
				s_solver_path1 += '/../bin/lindbladmpo.exe'
			else:
				s_cygwin_path1 = ''
				s_solver_path1 += '/../bin/lindbladmpo'
			if s_cygwin_path is None:
				s_cygwin_path = s_cygwin_path1
			if s_solver_path is None:
				s_solver_path = s_solver_path1
		return s_cygwin_path, s_solver_path

	@staticmethod
	def build_input_file(parameters: Dict) -> (str, str, str, str):
		""" Writing the input parameters from the input dictionary to txt file
		Args:
			parameters (dictionary): the arguments for the simulator
		Returns:
			(s_input_file, s_output_path, s_id_suffix): File name of solver input file, file prefix for solver output,
			 based on the user settings in the parameters dictionary (or default values if not assigned), and possibly
			appended with a unique id (if requested), the unique id string (possibly empty), and the file
			names suffix denoting the number of qubit, e.g. ".N=8".
		"""
		check_output = LindbladMPOSolver._check_argument_correctness(parameters)
		# check if there is a problem with the input, if "" returned there is no problem
		if check_output != "":
			raise Exception(check_output)
		# check if the user defined an input file name
		s_output_path = parameters.get("output_files_prefix", "")
		if s_output_path == "" or s_output_path[-1] in ('/', '.', '\\'):
			s_output_path += "lindblad"
		b_uuid = parameters.get("b_unique_id", False)
		s_id_suffix = ''
		if b_uuid:
			s_id_suffix = uuid.uuid4().hex
			print("Generating a unique id for this simulation: " + s_id_suffix)
			s_id_suffix = '.' + s_id_suffix
			s_output_path += s_id_suffix
		n_qubits = parameters.get("N")
		s_N_suffix = f".N={n_qubits}"
		s_input_file = s_output_path + s_N_suffix + ".input.txt"

		print("Creating solver input file:")
		print(s_input_file)
		AB_indices = False
		A_bond_indices = []
		B_bond_indices = []
		interactions = []
		if 'J' in parameters.keys():
			if type(parameters['J']) == np.ndarray:
				interactions.append('J')
		if 'J_z' in parameters.keys():
			if type(parameters['J_z']) == np.ndarray:
				interactions.append('J_z')
		file = open(s_input_file, "w")
		if len(interactions) == 2:
			if parameters['J'].shape == parameters['J_z'].shape:
				AB_indices = True
				for i in range(parameters['J'].shape[0]):
					for j in range(parameters['J'].shape[1]):
						if parameters['J'][i, j] != 0 or parameters['J_z'][i, j] != 0:
							A_bond_indices.append(i + 1)
							B_bond_indices.append(j + 1)
			else:
				raise Exception("J and J_z are not of the same size, aborting program.")
		elif len(interactions) == 1:
			AB_indices = True
			for i in range(parameters[interactions[0]].shape[0]):
				for j in range(parameters[interactions[0]].shape[1]):
					if parameters[interactions[0]][i, j] != 0:
						A_bond_indices.append(i + 1)
						B_bond_indices.append(j + 1)
		for key in parameters.keys():
			if key == "b_unique_id":
				pass
			elif key == "output_files_prefix":
				file.write(key + ' = ' + s_output_path + "\n")
			elif (key == 'J' or key == 'J_z') and type(parameters[key]) == np.ndarray and len(parameters[key]) > 1:
				# check if to create bond indices arrays
				not_first_value = False
				file.write(key + " = ")
				for i in range(len(A_bond_indices)):
					if not_first_value:
						file.write(",")
					file.write(str(parameters[key][A_bond_indices[i] - 1, B_bond_indices[i] - 1]))
					not_first_value = True
				file.write("\n")
			elif type(parameters[key]) == np.ndarray:
				file.write(key + " = ")
				for i in range(parameters[key].shape[0]):
					file.write(str(parameters[key][i]))
					if i + 1 != parameters[key].shape[0]:
						file.write(",")
				file.write("\n")
			else:
				file.write(key + " = " + str(parameters[key]).strip("[]") + "\n")
		if AB_indices:
			file.write("A_bond_indices = " + str(A_bond_indices).strip("[]") + "\n")
			file.write("B_bond_indices = " + str(B_bond_indices).strip("[]") + "\n")
		s_input_file = os.path.abspath(file.name).replace("\\","/")
		file.close()
		return s_input_file, s_output_path + s_N_suffix, s_id_suffix, s_N_suffix

	@staticmethod
	def execute(s_cygwin_path = None, s_solver_path = None, s_input_file = "") -> int:
		""" Execute the C++ solver
		Args:
			s_cygwin_path : the path of the cygwin bash terminal execution
			s_solver_path : the path of the simulator executable
			s_input_file : input file path, if empty then the simulator is run using default values.

		Returns:
			exit code : the exit code of the solver.
		"""
		s_cygwin_path, s_solver_path = LindbladMPOSolver.process_default_paths(s_cygwin_path, s_solver_path)
		if s_cygwin_path:
			call_string = s_cygwin_path + " --login -c \""
		else:
			call_string = ''
		call_string += s_solver_path
		if s_input_file:
			call_string += " input_file " + str(s_input_file)
		print("Executing solver with command:")
		print("\t" + call_string + "\n")

		process = subprocess.Popen(call_string, shell = True)
		exit_code = process.wait()
		print(f"Solver process terminated with exit code {exit_code}")
		return exit_code

	@staticmethod
	def load_output(s_output_path: str):
		result = {'1q': LindbladMPOSolver._read_1q_output_into_dict(s_output_path),
				  '2q': LindbladMPOSolver._read_2q_output_into_dict(s_output_path)}
		return result

	@staticmethod
	def _read_1q_output_into_dict(s_output_path: str) -> Dict:
		""" Reads the 1 qubit output file and returns a dictionary with the values in the format:
		key (tuple) = (qubit (int), Axis (string), time (float)) : value = expectation value (float)
		Args:
			s_output_path : file location
		Returns:
			result : the arguments for the simulator
		"""
		full_filename = s_output_path + ".obs.1q.dat"
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
	def _read_2q_output_into_dict(s_output_path: str) -> Dict:
		""" Reads the 2 qubit output file and returns a dictionary with the values in the format:
		key = (qubit1 (int), qubit2 (int), Axis1 (string), Axis2 (string) ,time (float)), value = expectation value (float)
		Args:
			s_output_path (string): file location
		Returns:
			result (dictionary): the arguments for the simulator
		"""
		full_filename = s_output_path + ".obs.2q.dat"
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
	def _get_number_of_qubits(dict_in: Dict) -> int:
		l_x_is_0 = False
		if "l_x" in dict_in:
			if dict_in["l_x"] == 0:
				l_x_is_0 = True
		if ("l_x" in dict_in) and ("l_y" in dict_in) and (not l_x_is_0):
			if LindbladMPOSolver._is_int(dict_in["l_x"]) and LindbladMPOSolver._is_int(dict_in["l_y"]):
				return dict_in["l_x"] * dict_in["l_y"]
		elif "N" in dict_in:
			if LindbladMPOSolver._is_int(dict_in["N"]):
				return dict_in["N"]
		return -1

	@staticmethod
	def _check_argument_correctness(dict_in: Dict) -> str:
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
		if ("N" not in dict_in) or ("t_final" not in dict_in) or ("tau" not in dict_in):
			check_msg += "Error 110: N, t_final and tau must be defined as they don't have default values\n"
			return check_msg
		for key in dict.keys(dict_in):
			if isinstance(dict_in[key], str) and "" == dict_in[key]:  # ignore empty entrances/space holders <"">
				continue
			flag_continue = False

			if key == "N":
				if not LindbladMPOSolver._is_int(dict_in[key]):
					check_msg += "Error 120: " + key + " should be an integer\n"
					continue
				if dict_in[key] <= 0:
					check_msg += "Error 130: " + key + " should be bigger/equal to 1 (integer)\n"
					continue

			elif (key == "t_final") or key == "tau":
				if not LindbladMPOSolver._is_float(dict_in[key]):
					check_msg += "Error 140: " + key + " is not a float\n"
					continue
				if dict_in[key] <= 0:
					check_msg += "Error 150: " + key + " must be bigger then 0\n"
					continue

			elif (key == "l_x") or (key == "l_y"):
				if not LindbladMPOSolver._is_int(dict_in[key]):
					check_msg += "Error 160: " + key + " should be an integer\n"
					continue
				if dict_in[key] < 0:
					check_msg += "Error 170: " + key + " should be bigger/equal to 1 (integer)\n"
					continue

			elif key == "output_step":
				if not LindbladMPOSolver._is_int(dict_in[key]):
					check_msg += "Error 180: " + key + " should be an integer\n"
					continue
				if dict_in[key] < 0:
					check_msg += "Error 190: " + key + " should be bigger/equal to 0 (integer)\n"
					continue

			elif (key == "h_x") or (key == "h_y") or (key == "h_z") or \
					(key == "g_0") or (key == "g_1") or (key == "g_2"):
				if LindbladMPOSolver._is_float(dict_in[key]):
					continue
				number_of_qubits = LindbladMPOSolver._get_number_of_qubits(dict_in)
				if number_of_qubits == -1:
					check_msg += "Error 200: " + key + " could not be validated because 'N' (or alternatively l_x, " \
													   "l_y) are not defined properly\n "
					continue
				if isinstance(dict_in[key], list):
					if len(dict_in[key]) != number_of_qubits:
						check_msg += "Error 210: " + key + " is not a float / N size list / numpy array (of floats)\n"
						continue
					for element in dict_in[key]:
						if not LindbladMPOSolver._is_float(element):
							check_msg += "Error 220: " + key + "is not a float / N size list / numpy array (of " \
															   "floats)\n "
							flag_continue = True
							break
					if flag_continue:
						continue
				elif isinstance(dict_in[key], np.ndarray):
					if (str((dict_in[key]).dtype).find("int") == -1) and (
							str((dict_in[key]).dtype).find("float") == -1):
						check_msg += "Error 230: " + key + " is not a float / N size list / numpy array (of floats)\n"
						continue
					if dict_in[key].size == 1:
						continue
					if (dict_in[key].shape[0] != number_of_qubits) or (dict_in[key].shape[0] != dict_in[key].size):
						check_msg += "Error 240: " + key + " is not a float / N size list / numpy array (of floats)\n"
						continue
				else:
					check_msg += "Error 250: " + key + " is not a float / N size list / numpy array (of floats)\n"
					continue

			elif (key == "J_z") or (key == "J"):
				if LindbladMPOSolver._is_float(dict_in[key]):
					continue
				number_of_qubits = LindbladMPOSolver._get_number_of_qubits(dict_in)
				if number_of_qubits == -1:
					check_msg += "Error 260: " + key + " could not be validated because 'N' (or alternatively l_x, " \
													   "l_y) are not defined properly\n"
					continue
				if isinstance(dict_in[key], list):
					if len(dict_in[key]) != number_of_qubits:
						check_msg += "Error 270: " + key + " should be a constant, or a square matrix (nested " \
														   "list/np.array) in the size of number_of_qubits^2 of floats\n"
						continue
					for lst in dict_in[key]:
						if not isinstance(lst, list):
							check_msg += "Error 280: " + key + "should be a constant, or a square matrix (nested " \
															   "list/np.array) in the size of number_of_qubits^2 of " \
															   "floats\n "
							flag_continue = True
							break
						if len(lst) != number_of_qubits:
							check_msg += "Error 290: " + key + "should be a constant, or a square matrix (nested " \
															   "list/np.array) in the size of number_of_qubits^2 of " \
															   "floats\n "
							flag_continue = True
							break
						for val in lst:
							if not LindbladMPOSolver._is_float(val):
								check_msg += "Error 300: " + key + "should be a constant, or a square matrix (nested " \
																   "list/np.array) in the size of number_of_qubits^2 " \
																   "of floats\n "
								flag_continue = True
								break
						if flag_continue:
							break
					if flag_continue:
						continue
				elif isinstance(dict_in[key], np.ndarray):
					if (str((dict_in[key]).dtype).find("int") == -1) and (
							str((dict_in[key]).dtype).find("float") == -1):
						check_msg += "Error 310: " + key + "should be a constant, or a square matrix (nested " \
														   "list/np.array) in the size of number_of_qubits^2 of " \
														   "floats\n "
						continue
					if dict_in[key].size == 1:
						continue
					if dict_in[key].shape[0] != number_of_qubits:
						check_msg += "Error 320: " + key + "should be a constant, or a square matrix (nested " \
														   "list/np.array) in the size of number_of_qubits^2 of " \
														   "floats\n "
						continue
					if dict_in[key].shape[0] ** 2 != dict_in[key].size:
						check_msg += "Error 330: " + key + "should be a constant, or a square matrix (nested " \
														   "list/np.array) in the size of number_of_qubits^2 of " \
														   "floats\n "
						continue
				else:
					check_msg += "Error 340: " + key + " should be a constant, or a square matrix (nested " \
													   "list/np.array) in the size of number_of_qubits^2 of floats\n"
					continue

			elif key == "init_Pauli_state":
				if not isinstance(dict_in[key], str):
					check_msg += "Error 350: " + key + " is not a string\n"
					continue
				str_as_lst = list(dict_in[key])
				if len(str_as_lst) != 2:
					check_msg += "Error 360: " + key + " is a string but not of length 2 !\n"
					continue
				if (str_as_lst[0] != '+') and (str_as_lst[0] != "-"):
					check_msg += "Error 370: " + key + " can only be one of these: +x, -x, +y, -y, +z, -z\n"
					continue
				if (str_as_lst[1] != "x") and (str_as_lst[1] != "y") and (str_as_lst[1] != "z"):
					check_msg += "Error 380: " + key + " can only be one of these: +x, -x, +y, -y, +z, -z\n"
					continue

			elif ((key == "b_periodic_x") or (key == "b_periodic_y") or (key == "b_force_rho_trace") or (
					key == "b_force_rho_hermitian") or (key == "b_unique_id") or (key == "b_save_final_state")):
				if not isinstance(dict_in[key], bool):
					check_msg += "Error 390: " + key + " should be a boolean True or False as its a switch\n"
					continue

			elif key == "trotter_order":
				if not LindbladMPOSolver._is_int(dict_in[key]):
					check_msg += "Error 400: " + key + " should be 2, 3 or 4\n"
					continue
				if (dict_in[key] != 2) and (dict_in[key] != 3) and (dict_in[key] != 4):
					check_msg += "Error 401: " + key + " should be 2, 3 or 4\n"
					continue

			elif (key == "max_dim") or (key == "max_dim_rho"):  # int
				if not LindbladMPOSolver._is_int(dict_in[key]):
					check_msg += "Error 410: " + key + " should be an integer\n"
					continue

			elif (key == "cut_off") or (key == "cut_off_rho"):
				if not LindbladMPOSolver._is_float(dict_in[key]):
					check_msg += "Error 420: " + key + " is not a small float format (ae-b where a and b are numbers)\n"
					continue

			elif key == "output_files_prefix":
				pass
			elif key == "1q_components":
				x_c = 0
				y_c = 0
				z_c = 0
				if not isinstance(dict_in[key], list):
					check_msg += "Error 430: " + key + " should get a list of sizes 1,2,3 with x,y,z)\n"
					continue
				if len(dict_in[key]) > 3:
					check_msg += "Error 440: " + key + " should get a list of sizes 1,2,3 with x,y,z)\n"
					continue
				for val in dict_in[key]:
					if not isinstance(val, str):
						check_msg += "Error 441: " + key + " only gets x,y,z (or a subset)\n"
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
						check_msg += "Error 450: " + key + " only gets x,y,z (or a subset)\n"
						flag_continue = True
						break
				if flag_continue:
					continue
				if (x_c > 1) or (y_c > 1) or (z_c > 1):
					check_msg += "Error 460: " + key + " only gets x,y,z (or a subset)\n"
					continue

			elif key == "1q_indices":
				if dict_in[key] != "":
					if not isinstance(dict_in[key], list):
						check_msg += "Error 470: " + key + " should be an integer list (1,2,3,4..)\n"
						continue
					number_of_qubits = LindbladMPOSolver._get_number_of_qubits(dict_in)
					if number_of_qubits == -1:
						check_msg += "Error 480: " + key + "could not be validated because 'N' (or alternatively l_x," \
														   " l_y) are not defined properly\n "
						continue
					for element in dict_in[key]:
						if not LindbladMPOSolver._is_int(element):
							check_msg += "Error 490: " + key + " should be an integer list (1,2,3,4..)\n"
							flag_continue = True
							break
						if element >= number_of_qubits:
							check_msg += "Error 500: " + key + " should be an integer list listing qubits, therefore " \
															   "number range is 0 to 'num of qubits'-1\n"
							flag_continue = True
							break
					if flag_continue:
						continue
					if len(dict_in[key]) > number_of_qubits:
						check_msg += "Error 510: " + key + " 's length should be smaller/equal then the amount of " \
														   "qubits\n "
						continue
					if not len(set(dict_in[key])) == len(dict_in[key]):
						check_msg += "Error 520: " + key + " 's List does not contains all unique elements"
						continue

			elif key == "2q_components":
				if not isinstance(dict_in[key], list):
					check_msg += "Error 530: " + key + "only receives xx,yy,zz,xy,xz,yz (or a subset) as a strings " \
													   "list\n"
					continue
				if len(dict_in[key]) > 6:
					check_msg += "Error 540: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
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
						check_msg += "Error 550: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
						flag_continue = True
						break
				if flag_continue:
					continue
				for check_val in check_me:
					if check_val > 1:
						check_msg += "Error 560: " + key + " only receives xx,yy,zz,xy,xz,yz (or a subset)\n"
						flag_continue = True
						break
				if flag_continue:
					continue

			elif key == "2q_indices":  # expecting an integer tuples list
				if not isinstance(dict_in[key], list):
					check_msg += "Error 570: " + key + " should be an list of tuples of size 2, containing integer\n"
					continue
				number_of_qubits = LindbladMPOSolver._get_number_of_qubits(dict_in)
				if number_of_qubits == -1:
					check_msg += "Error 580: " + key + " could not be validated because 'N' (or alternatively l_x, " \
													   "l_y) are not defined properly\n"
					continue
				for tup in dict_in[key]:
					if not isinstance(tup, tuple):
						check_msg += "Error 590: " + key + " should be an list of tuples of size 2, containing " \
														   "integer\n "
						flag_continue = True
						break
					if ((not LindbladMPOSolver._is_int(tup[0])) or (not LindbladMPOSolver._is_int(tup[1])) or (
							len(tup) != 2)):
						check_msg += "Error 600: " + key + " should be an list of tuples of size 2, containing " \
														   "integers\n "
						flag_continue = True
						break
					if (tup[0] >= number_of_qubits) or (tup[1] >= number_of_qubits):
						check_msg += "Error 610: " + key + " should be an list of tuples of size 2, containing " \
														   "integers smaller/equal then the total number of qubits\n "
						flag_continue = True
						break
				if flag_continue:
					continue

				if len(dict_in[key]) > number_of_qubits ** 2:
					check_msg += "Error 620: " + key + " 's length should be smaller then the amount of qubits^2\n"
					continue
				if not len(set(dict_in[key])) == len(dict_in[key]):
					check_msg += "Error 630: " + key + " 's List does not contains all unique elements"
					continue

			else:
				check_msg += "Error: invalid key (parameter) inserted: " + key + "\n"
		# End of: "for key in dict.keys(dict_in)"
		# More cross-parameter checks:
		if ("t_final" in dict_in) and ("tau" in dict_in):
			if (LindbladMPOSolver._is_float(dict_in["tau"])) and (LindbladMPOSolver._is_float(dict_in["t_final"])):
				if (dict_in["tau"] > 0) and (dict_in["t_final"] > 0):
					if dict_in["tau"] > dict_in["t_final"]:
						check_msg += "Error 640: t_final (total time) is smaller then tau (time step for time " \
									 "evolution)\n "
					elif "output_step" in dict_in:
						if LindbladMPOSolver._is_int(dict_in["output_step"]):
							if dict_in["output_step"] > 0:
								if dict_in["output_step"] * dict_in["tau"] > dict_in["t_final"]:
									check_msg += "Error 650: Output_step multiplied by tau, is bigger then t_final (" \
												 "output_step in units of tau, times tau is bigger then the " \
												 "simulation time)\n "
		return check_msg
