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

# Simulation parameters
n_qubits = 8
b_save_figures = True
fontsize = 22

h_x = 0. * np.random.randn(n_qubits)
h_x[int(n_qubits / 2)] = .5
h_y = 0. * np.random.randn(n_qubits)

h_z = 0. * np.random.randn(n_qubits)
g_0 = .1 * np.random.rand(n_qubits)
J = 1
t_final = 5
tau = .05

# Create the parameters dictionary
s_file_prefix = f"disordered chain"
solver_params = {'tau': tau, 't_final': t_final, 'max_dim_rho': 80, 'N': n_qubits, 'b_unique_id': False,
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
plot_space_time(data, t_steps, t_ticks, qubits, fontsize = fontsize, s_file_prefix = s_output_path + s_file_prefix)
plt.show()
