# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
from driving_routines import *

"""
An example of a large scale research project for the dynamics qubit with XY (flip-flop) interactions
with different topologies (connectivity), in which an edge qubit is being driven on resonance
in the lab frame, and the qubits have alternating frequencies between nearest neighbors.
The problem is simulated in the frame rotating with the driven qubit's frequency (after a rotating
wave approximation), with the qubits initialized along the z axis of the Bloch sphere
(the pure product ground state), and subject to energy relaxation and decoherence.     
"""

# Plotting parameters
b_save_figures = True		# Whether to save the plotted figures to disk
b_save_to_db = True			# Whether to save to a .csv file an entry describing the simulation
fontsize = 22 				# Font size of the figures

# Simulation parameters to vary - all of these are the metadata (database) entries
N = 9 						# Number of qubits. Note that the topologies are implemented only
# for some qubit numbers - see in the following.
topology = 'chain.E'  		# Options: 'chain.M', 'chain.E', 'ring', 'plaquette'
# Setting 'chain.M', the system is a linear 1D chain, with the Middle qubit driven by an x term,
# Setting 'chain.E', the system is a linear 1D chain, with the Edge (0) qubit driven by an x term,
# Setting 'ring', the system is a closed ring, with qubit 0 driven by an x term,
# Setting 'plaquette', the system is a closed ring with two additional edge qubits at opposite side,
# 	with one edge qubit (0) driven by an x term.
solver = 'scipy'  			# Options: 'mpo', 'scipy'
# Setting 'mpo' the MPO executable is invoked for the solution, while 'scipy' uses qiskit-dynamics
# and its scipy solver

load_unique_id = ''			# uuid of a previous simulation, in order to start the current simulation
# with the initial condition being the final state of the loaded simulation
t_init = 0.					# Initial time, can be set to the final time of the loaded simulation,
# which is convenient for plotting but not essential for the solution
t_final = 1. 				# Final time of ths simulation
tau = .02 					# Discrete time step of the Trotterization. Should be much smaller than
# all periods determined by the Hamiltonian parameters below.

# Model parameters
J_amp = 1.					# XY coupling strength * 2 * pi
h_x_amp = 1. 				# x-driving amplitude of the 0 qubit * 2 * pi
h_z_amp = 0. 				# On-site z terms amplitude according to the alternating pattern
g_0_amp = .1 				# Rate of energy relaxation towards the ground state
g_2_amp = .0 				# Rate of decoherence (dephasing) terms in the Lindbladian

if topology == 'plaquette' or topology == 'ring':
	mpo_mapping = 'A'  		# Options: 'A', 'B' (with 'B' only implemented for some qubit numbers)
else:
	mpo_mapping = ''
# The mpo_mapping option indicates the ordering of qubits in the topologies. It makes a difference
# in the MPO solution only, because of the data structure holding the solution ansatz.

# See below for further important numerical accuracy parameters according to the solver type.

# Create the metadata (database) dictionary, which defines our simulation. The default one provides
# empty values for all necessary fields
sim_metadata = DEF_METADATA.copy()
sim_metadata.update({
	'topology': topology, 'N': N, 'solver': solver, 't_init': t_init, 't_final': t_final, 'tau': tau,
	'J_amp': J_amp, 'h_x_amp': h_x_amp, 'h_z_amp': h_z_amp, 'g_0_amp': g_0_amp, 'g_2_amp': g_2_amp,
	'load_unique_id': load_unique_id})

# Now define the solver-specific parameters, and update the metadata
if solver == 'mpo':
	max_dim_rho = 60  				# Maximal bond dimension
	cut_off_rho = 1e-14  			# Cutoff parameter of the Schmidt decomposition singular values
	force_rho_Hermitian_step = 5  	# Every how many steps to force the density matrix to be Hermitian

	sim_metadata.update({'mpo_mapping': mpo_mapping,
						 'cut_off_rho': cut_off_rho,
						 'max_dim_rho': max_dim_rho,
						 'force_rho_Hermitian_step': force_rho_Hermitian_step})
elif solver == 'scipy':
	method = 'RK45'  				# Scipy solver, for example: 'RK45', 'DOP853'
	atol = 1e-8  					# Absolute tolerance parameter
	rtol = 1e-8  					# Relative tolerance parameter
	sim_metadata.update({'mpo_mapping': mpo_mapping, 'method': method, 'atol': atol, 'rtol': rtol})

# Use the metadata dictionary to build parameters dictionary, create a solver, save database entry,
# solve, save the output results, and plot simulation figures.
solve_simulation(sim_metadata, fontsize, b_save_to_db = b_save_to_db, b_save_figures = b_save_figures)
