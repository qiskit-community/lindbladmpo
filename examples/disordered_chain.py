# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from lindbladmpo.LindbladMPOSolver import *
from lindbladmpo.plot_routines import *

s_output_path = os.path.abspath('./output') + '/'
if not os.path.exists(s_output_path):
	os.mkdir(s_output_path)
s_file_prefix = "chain"

# Simulation parameters
rand_seed = 1
n_qubits = 9
b_unique_id = True
b_save_figures = False
fontsize = 22
h_z_half_width = 0.

np.random.seed(rand_seed)
h_x = np.zeros(n_qubits, float)
h_x[int(n_qubits / 2)] = 5.

h_z = -h_z_half_width + 2 * h_z_half_width * np.random.rand(n_qubits)  # np.ones(n_qubits)  # np.random.rand(n_qubits)  #
g_0 = .1 * np.ones(n_qubits)
J = 1
t_final = 10
tau = .01

# Create the parameters dictionary
solver_params = {'tau': tau, 't_final': t_final, 'cut_off_rho': 1e-16,
				 'max_dim_rho': 100, 'N': n_qubits, 'b_unique_id': b_unique_id,
				 'h_x': h_x, 'h_z': h_z, 'g_0': g_0, 'J': J, 'l_x': n_qubits, 'l_y': 1,
				 'output_files_prefix': s_output_path + s_file_prefix}
# Initialize class arguments - the parameters, cygwin path, and MPO executable path
s_cygwin_path = None
s_solver_path = None
solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
# Execute simulator
solver.solve()
# Prepare plot data
data, t_steps, t_ticks, qubits = prepare_plot_data(solver)
# And plot one figure. By default it will be saved to a file.
plot_space_time(data, t_steps, t_ticks, qubits, fontsize = fontsize, s_file_prefix = s_output_path + s_file_prefix,
				b_save_figures = b_save_figures)
plt.show()
