# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Routines for managing the research project running multiple simulations.
"""

import os.path
from topologies import *
from output_routines import *
from lindbladmpo.plot_routines import *
from lindbladmpo.examples.simulation_building.LindbladMatrixSolver import LindbladMatrixSolver


DEF_METADATA = {'unique_id': '', 'topology': '', 'N': '', 'solver': '',
				't_init': '', 't_final': '', 'tau': '',
				'J_amp': '', 'h_x_amp': '', 'h_z_amp': '', 'g_0_amp': '', 'g_2_amp': '',
				'unused1': '', 'max_dim_rho': '', 'cut_off_rho': '',
				'mpo_mapping': '', 'method': '', 'atol': '', 'rtol': '',
				'load_unique_id': '', 'force_rho_Hermitian_step': '', 'unused2': ''}
"""Default metadata and db parameters for the simulations in this project."""

S_DB_FILENAME = 'simulations.csv'
"""File name to use for the database of the project simulations."""

S_FILE_PREFIX = 'sim'
"""Prefix for the file names of all solver output files and plots of the project simulations."""


def solve_simulation(sim_metadata: dict, fontsize = 20, b_save_to_db = True, b_save_figures = True):
	topology = sim_metadata['topology']
	N = sim_metadata['N']
	s_solver = sim_metadata['solver']
	t_init = sim_metadata['t_init']
	t_final = sim_metadata['t_final']
	tau = sim_metadata['tau']
	J_amp = sim_metadata['J_amp']
	h_x_amp = sim_metadata['h_x_amp']
	h_z_amp = sim_metadata['h_z_amp']
	g_0_amp = sim_metadata['g_0_amp']
	g_2_amp = sim_metadata['g_2_amp']
	mpo_mapping = sim_metadata['mpo_mapping']
	load_unique_id = sim_metadata['load_unique_id']
	force_rho_Hermitian_step = sim_metadata['force_rho_Hermitian_step']

	# -------------------------------------------------------------------------
	# The following parameters are fixed for the project and kept hard-coded.
	_1q_components = ['x', 'y', 'z']
	_2q_components = ['xy']
	_1q_plot_components = ['x', 'z']
	_2q_plot_components = ['xy']

	# _1q_components = ['x', 'y', 'z']
	# _2q_components = ['zz', 'zx', 'xx', 'xy', 'zy', 'yy']
	# _1q_plot_components = ['x', 'y', 'z']
	# _2q_plot_components = ['zz', 'xx', 'xy']

	init_pauli_state = []
	_1q_indices = []; _2q_indices = []
	r_qubits = range(N)
	for i in r_qubits:
		if load_unique_id == '':
			init_pauli_state.append('+z')  # Initialize all qubits to the ground state
		_1q_indices.append(i)  # Request 1Q observables all qubits
	if topology == 'chain.M':
		# The qubits to plot at the end of the simulation are by default not all qubits.
		if N < 3:
			raise Exception("At least 3 qubits must be in the driven chain.")
		i_driven_q = int(N / 2)
		_1q_plot_indices = [i_driven_q, i_driven_q - 1, i_driven_q + 1]
		_2q_plot_indices = [(0, i_driven_q)]
		if N >= 5:
			_1q_plot_indices.extend([0, N - 1])
			_2q_plot_indices.extend([(1, i_driven_q), (i_driven_q, N - 2)])
		_2q_plot_indices.append((i_driven_q, N - 1))
	elif topology == 'chain.E':
		# The qubits to plot at the end of the simulation are by default not all qubits.
		if N < 3:
			raise Exception("At least 3 qubits must be in the driven chain.")
		i_driven_q = 0
		_1q_plot_indices = [0, 1, 2]
		if N >= 6:
			_1q_plot_indices.append((N - 3))
		if N >= 5:
			_1q_plot_indices.append((N - 2))
		if N >= 4:
			_1q_plot_indices.append((N - 1))
		_2q_plot_indices = [(0, 1)]
		if N >= 4:
			_2q_plot_indices.append((0, 2))
		if N >= 6:
			_2q_plot_indices.append((0, 3))
		if N >= 7:
			_2q_plot_indices.append((0, N - 3))
		if N >= 5:
			_2q_plot_indices.append((0, N - 2))
		_2q_plot_indices.append((0, N - 1))
	elif topology == 'plaquette' or topology == 'ring':
		if N < 6:
			raise Exception("At least 6 qubits must be in the driven plaquette or ring.")
		i_driven_q = 0
		if mpo_mapping == 'A':
			_1q_plot_indices = [0, 2, 3]
			_2q_plot_indices = [(0, 1), (0, 2), (0, 3)]
			if N >= 8:
				_1q_plot_indices.extend([N - 4, N - 3])
				_2q_plot_indices.extend([(0, N - 4), (0, N - 3)])
			_1q_plot_indices.append(N - 1)
			_2q_plot_indices.append((0, N - 1))
		elif mpo_mapping == 'B':
			if topology == 'ring':
				_1q_plot_indices = [0, 1, N - 1]
				_2q_plot_indices = [(0, 1), (0, 2), (0, N - 1)]
				left_Q = int(N / 2) - 2
			else:  # plaquette
				_1q_plot_indices = [0, 2, N - 1]
				_2q_plot_indices = [(0, 1), (0, 2), (0, N - 1)]
				left_Q = int(N / 2) - 1
			if N >= 8:
				_1q_plot_indices.extend([left_Q, left_Q + 3])
				_2q_plot_indices.extend([(0, left_Q), (0, left_Q + 3)])
			_1q_plot_indices.append(left_Q + 2)
			_2q_plot_indices.append((0, left_Q + 2))
		else:
			raise Exception(f"Unknown mpo_mapping '{mpo_mapping}'.")
	else:
		raise Exception(f"Unknown topology '{topology}'.")
	if len(_2q_components):
		# Add 2Q observables for driven qubit with all others. For consistency fix the first index
		# to be the lower of the two.
		for i in range(0, N):
			for j in range(0, N):
				if i != j:
					_2q_indices.append((i, j))

	s_output_path = os.path.abspath('./output') + '/'  # use the absolute path of the current file
	s_data_path, s_plot_path = generate_paths(s_output_path)
	s_data_path += S_FILE_PREFIX
	load_files_prefix = ''
	if load_unique_id != '':
		load_files_prefix = s_data_path + '.' + load_unique_id

	s_coupling_map = f"{N}.{topology}"
	if 'chain' not in topology:
		s_coupling_map += f".{mpo_mapping}"

	s_plot_path += S_FILE_PREFIX + f".N={s_coupling_map}"

	# -------------------------------------------------------------------------
	# Execution section

	_2_pi = 2 * np.pi
	h_z_pattern = h_z_patterns[s_coupling_map]
	coupling_map = coupling_maps[s_coupling_map]
	h_x = np.zeros(N)
	h_x[i_driven_q] = h_x_amp * _2_pi
	h_z = h_z_amp * _2_pi * np.asarray(h_z_pattern)
	g_0 = g_0_amp * np.ones(N)
	g_2 = g_2_amp * np.ones(N)
	J = J_amp * _2_pi
	J_ij = np.zeros((N, N))
	for bond in coupling_map:
		J_ij[bond[0], bond[1]] = J

	plot_topology(N, topology, s_coupling_map, b_save_figures, s_plot_path,
				  'chain' not in topology, h_z_amp != 0)

	b_save_final_state = True
	if N > 10 and s_solver == 'scipy':
		b_save_final_state = False
	solver_params = {'N': N, 'b_unique_id': True, 'metadata': str(sim_metadata),
					 't_init': t_init, 't_final': t_final, 'tau': tau,
					 'h_x': h_x, 'h_z': h_z, 'g_0': g_0, 'g_2': g_2, 'J': J_ij,
					 'init_pauli_state': init_pauli_state, 'load_files_prefix': load_files_prefix,
					 '1q_components': _1q_components, '1q_indices': _1q_indices,
					 '2q_components': _2q_components, '2q_indices': _2q_indices,
					 'output_files_prefix': s_data_path,
					 'b_save_final_state': b_save_final_state,
					 'output_step': 1, 'force_rho_hermitian_step': force_rho_Hermitian_step}

	# Create the solver, and update solver-specific parameters
	if s_solver == 'mpo':
		solver = LindbladMPOSolver()
		solver_params.update({'max_dim_rho': sim_metadata['max_dim_rho'],
							  'cut_off_rho': sim_metadata['cut_off_rho']})
	elif s_solver == 'scipy':
		solver = LindbladMatrixSolver()
		solver_params.update({'method': sim_metadata['method'],
							  'atol': sim_metadata['atol'], 'rtol': sim_metadata['rtol']})
	else:
		raise Exception("Solver type is unsupported.")

	# Verify parameters and create the solver input file. This will create also a unique id fot the
	# solver output files, that we store in the db.
	solver.build(solver_params)
	sim_metadata['unique_id'] = solver.parameters['unique_id']

	if b_save_to_db:
		s_db_path = s_output_path + S_DB_FILENAME
		save_to_db(s_db_path, sim_metadata)

	# Solve, save final state and observable data files.
	solver.solve()

	# Plot figures
	result = solver.result
	parameters = solver.parameters
	s_file_prefix = s_plot_path + solver.s_id_suffix

	# Plot 1Q space-time diagrams, and a few select time curves
	for s_obs_name in _1q_plot_components:
		plot_full_1q_space_time(parameters, result, s_obs_name, fontsize = fontsize,
								b_save_figures = b_save_figures, s_file_prefix = s_file_prefix)
		plot_1q_obs_curves(parameters, result, s_obs_name, _1q_plot_indices, fontsize = fontsize,
						   b_save_figures = b_save_figures, s_file_prefix = s_file_prefix)

	# Plot 2Q connected correlation functions as space-time diagrams, and a few select time curves
	for s_obs_name in _2q_plot_components:
		plot_2q_obs_curves(parameters, result, s_obs_name, _2q_plot_indices, fontsize = fontsize,
						   b_save_figures = b_save_figures, s_file_prefix = s_file_prefix)
		plot_full_2q_correlation_matrix(solver.parameters, solver.result, s_obs_name,
										t_final, fontsize = fontsize,
										b_save_figures = b_save_figures, s_file_prefix = s_file_prefix)

	# Plot the second Renyi entropy (minus the log of the purity loss).
	plot_global_obs_curve(parameters, result, 'S_2', fontsize = 16,
						  b_save_figures = b_save_figures, s_file_prefix = s_file_prefix)
	plt.show()
