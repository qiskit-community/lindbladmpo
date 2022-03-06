# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from lindbladmpo.LindbladMPOSolver import *
from lindbladmpo.plot_routines import *

"""
An example of the dynamics of a 1D qubit chain with XY (flip-flop) interactions, in which the
center qubit is being driven on resonance in the lab frame, and the qubits have random frequencies.
The problem is simulated in the frame rotating with the driven qubit's frequency (after a rotating
wave approximation), with the qubits initialized along the z axis of the Bloch sphere
(the pure product ground state), and subject to energy relaxation and decoherence.      
"""

s_output_path = os.path.abspath('./output') + '/'  # All solver output files will be written here
if not os.path.exists(s_output_path):
	os.mkdir(s_output_path)
s_file_prefix = "chain"  # All solver output files will begin with this prefix

# Input and output parameters
rand_seed = 1  			# Seed for the numpy random function
b_unique_id = True  	# Whether to generate a uuid for the file names of the solver output
b_save_figures = True  	# Whether to save .png files plotted at the end of simulation
fontsize = 22  			# Font size in the figures

# Primary model parameters - vary the values of those and compare the results.

n_qubits = 21 			# Number of qubits to simulate. Explore different numbers of qubits
# (with tens and even hundreds). Use an odd number to make the setup reflection-symmetric,
# as the middle qubit has an x term below.
h_z_amp = 10.  			# Qubit frequencies will be random (uniformly distributed) in the interval
# 2 * pi * [-h_z_amp, h_z_amp].
# For h_z_amp the qubits ar resonant, and there is no localization of the driven qubit's excitations.
# For h_z_amp large enough, the excitations of the driven qubit will be inhibited from propagating.
t_final = 1.  			# Final simulation time
tau = 0.01  			# The discrete integration (Trotterization) time step size.
# It should be significantly smaller than the period of all coherent dynamics.
max_dim_rho = 60  		# This is the maximal bond dimension at each qubit.

# Other model parameters
g_0 = .1  				# Energy relaxation rate, uniform for all qubits
g_2 = .1  				# Decoherence (dephasing) rate, uniform for all qubits
J = 1. * 2 * np.pi  	# XY coupling coefficient, uniform for all qubit bonds.
h_x_amp = 0.5  			# The center qubit in the chain has an x-axis term in the Hamiltonian
# with this amplitude (* 2 * pi).

# Initialize the h_z and h_x parameter arrays
h_x = np.zeros(n_qubits, float)
h_x[int(n_qubits / 2)] = 0.5 * 2 * np.pi
np.random.seed(rand_seed)
h_z = -h_z_amp + 2 * h_z_amp * np.random.rand(n_qubits)

# Create the parameters dictionary
solver_params = {
	'tau': tau,
	't_final': t_final,
	'cut_off_rho': 1e-14,
	'max_dim_rho': max_dim_rho,
	'N': n_qubits,
	'b_unique_id': b_unique_id,
	'h_x': h_x,
	'h_z': h_z,
	'g_0': g_0,
	'g_2': g_2,
	'J': J,
	'1q_components': ['x', 'y', 'z'],
	'2q_components': [],
	'init_pauli_state': '+z',
	'output_files_prefix': s_output_path + s_file_prefix
}
# Initialize class arguments - the parameters, cygwin path, and MPO executable path
solver = LindbladMPOSolver(solver_params)
# Execute simulator
solver.solve()
# Plot one figure, a space-time diagram of the Z Pauli observable expectation value.
plot_full_1q_space_time(solver.parameters, solver.result, 'z', fontsize = fontsize,
						b_save_figures = b_save_figures, s_file_prefix = s_output_path + s_file_prefix)
plt.show()
