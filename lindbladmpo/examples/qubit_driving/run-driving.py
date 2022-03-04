# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
from driving_routines import *

# Plotting parameters
b_save_figures = True
b_save_to_db = True
fontsize = 22

# Simulation parameters to vary - all of these are the metadata (database) entries
N = 8
topology = 'plaquette'  # Options: 'chain.M', 'chain.E', 'plaquette', 'ring'
solver = 'scipy'  # Options: 'mpo', 'scipy'
if topology == 'plaquette' or topology == 'ring':
	mpo_mapping = 'A'  # Options: 'A', 'B'. The 'B' mapping is only implemented for select qubit numbers.
else:
	mpo_mapping = ''

J_amp = 1.
h_x_amp = 1.
h_z_amp = 0.
g_0_amp = .1
g_2_amp = .0
t_init = 0.
t_final = 1.
load_unique_id = ''
tau = .02

# Create the metadata (database) dictionary, which defines our simulation. The default one provides
# empty values for all necessary fields
sim_metadata = DEF_METADATA.copy()
sim_metadata.update({
	'topology': topology, 'N': N, 'solver': solver, 't_init': t_init, 't_final': t_final, 'tau': tau,
	'J_amp': J_amp, 'h_x_amp': h_x_amp, 'h_z_amp': h_z_amp, 'g_0_amp': g_0_amp, 'g_2_amp': g_2_amp,
	'load_unique_id': load_unique_id})

# Now define the solver-specific parameters, and update the metadata
if solver == 'mpo':
	max_dim_rho = 60  # Maximal bond dimension
	cut_off_rho = 1e-14  # Cutoff parameter of the Schmidt decomposition singular values
	force_rho_Hermitian_step = 5  # Every how many steps to force the density matrix to be Hermitian

	sim_metadata.update({'mpo_mapping': mpo_mapping,
						 'cut_off_rho': cut_off_rho,
						 'max_dim_rho': max_dim_rho,
						 'force_rho_Hermitian_step': force_rho_Hermitian_step})
elif solver == 'scipy':
	method = 'RK45'  # Scipy solver, for example: 'RK45', 'DOP853'
	atol = 1e-8  # Absolute tolerance parameter
	rtol = 1e-8  # Relative tolerance parameter
	sim_metadata.update({'mpo_mapping': mpo_mapping, 'method': method, 'atol': atol, 'rtol': rtol})

# Use the metadata dictionary to build parameters dictionary, create a solver, save database entry,
# solve, save the output results, and plot simulation figures.
solve_simulation(sim_metadata, fontsize, b_save_to_db = b_save_to_db, b_save_figures = b_save_figures)

tmp = 2  # Put a breakpoint here if desired
