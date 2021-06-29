# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from MPO_Lindblad_solver.MPOLindbladSolver import MPOLindbladSolver
from qiskit.visualization.gate_map import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from MPO_Lindblad_solver.temp_gate_map import *

# Simulation parameters
n_qubits = 8
b_plaquette = True
# s_path = "C:/Users/galvz/PycharmProjects/sim_func/"
# s_executable = "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator/lindblad.exe"
s_path = "C:/temp/"
s_executable = "/cygdrive/c/Users/HaggaiLanda/gitprojects/MPO-Lindblad-solver/src/lindblad.exe"
b_save_figures = True
fontsize = 22

h_x = 0. * np.random.randn(n_qubits)
h_x[int(n_qubits / 2)] = 5.
h_y = 0. * np.random.randn(n_qubits)

h_z = 0. * np.random.randn(n_qubits)
g_1 = (.1 * np.random.rand(n_qubits)).tolist()
J = 1
t_final = .5
tau = .05

mpl_data = {}
c_map = {}
if not b_plaquette:
	l_x = n_qubits
	l_y = 1
	J_param = J
else:
	c_map[27] = [[0, 1], [1, 0], [1, 2], [1, 4], [2, 1], [2, 3], [3, 2], [3, 5], [4, 1], [4, 7], [5, 3], [5, 8], [6, 7], [7, 4], [7, 6], [7, 10], [8, 5], [8, 9], [8, 11], [9, 8], [10, 7], [10, 12], [11, 8], [11, 14], [12, 10], [12, 13], [12, 15], [13, 12], [13, 14], [14, 11], [14, 13], [14, 16], [15, 12], [15, 18], [16, 14], [16, 19], [17, 18], [18, 15], [18, 17], [18, 21], [19, 16], [19, 20], [19, 22], [20, 19], [21, 18], [21, 23], [22, 19], [22, 25], [23, 21], [23, 24], [24, 23], [24, 25], [25, 22], [25, 24], [25, 26], [26, 25]]
	c_map[18] = [[0, 1], [1, 2], [1, 17], [2, 4], [4, 3], [4, 5], [5, 7], [6, 7], [7, 8], [8, 9], [9, 10], [9, 11],
			 [11, 12], [12, 13], [12, 14], [14, 15], [15, 16], [15, 17]]
	mpl_data[18] = [[0, 3], [1, 3], [1, 2], [1, 0], [1, 1],
					[2, 1], [3, 0], [3, 1], [3, 2], [3, 3],
					[4, 3], [3, 4], [3, 5], [3, 6], [2, 5],
					[1, 5], [1, 6], [1, 4]]

	c_map[8] = [[0, 1], [1, 2], [2, 3], [3, 5], [4, 5], [5, 6], [6, 7], [7, 1]]
	mpl_data[8] = [[0, 1], [1, 1], [1, 0], [2, 0], [3, 1],
					[2, 1], [2, 2], [1, 2]]

	J_ij = np.zeros((n_qubits, n_qubits))
	for ab in c_map[n_qubits]:
		J_ij[ab[0], ab[1]] = J
	J_param = J_ij
	l_x = 0
	l_y = 1

# create the parameters dictionary
solver_params = {'tau': tau, 't_final': t_final, 'max_dim_rho': 60, 'N': n_qubits, 'b_unique_id': False,
				 'h_x': h_x, 'h_z': h_z, 'g_1': g_1, 'J': J_param, 'l_x': l_x, 'l_y': l_y,
				 'input_file_prefix': s_path, 'output_file_prefix': s_path + "MPO"}
# initialize class - parameters, cygwin_path, simulator_path
solver = MPOLindbladSolver(solver_params, "C:/cygwin64/bin/bash.exe", s_executable)

# execute simulator
solver.solve()

n_t_steps = int(t_final/tau) + 1
data = np.zeros(shape=(n_qubits, n_t_steps), dtype = float)
i = 0
for key in solver.result['1q']:
	if key[1] == 'Z':
		data[key[0] - 1, i] = solver.result['1q'][key]
		if key[0] == n_qubits:
			i += 1
t_steps = np.arange(0, n_t_steps, int(n_t_steps / 10))
t_ticks = np.round(t_steps * tau, 5)
r_qubits = range(n_qubits)

if b_plaquette:
	plot_gate_map_raw(n_qubits, mpl_data, c_map[n_qubits], qubit_size = 28)

fig, ax = plt.subplots(figsize = (14, 9))
plt.rcParams.update({'font.size': fontsize})
im = ax.imshow(data, interpolation='none', aspect = 'auto')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
# plt.colorbar(im)
ax.set_xlabel('$t$', fontsize = fontsize)
ax.set_xticks(t_steps)
ax.set_xticklabels(t_ticks, fontsize = fontsize)
ax.set_yticks(r_qubits)
ax.set_yticklabels(r_qubits, fontsize = fontsize)
ax.set_ylabel('qubits', fontsize = fontsize)
ax.set_title('$\\langle\sigma_z(t)\\rangle$')
s_file = f"sigma_z(t), N = {n_qubits}, disordered 1D lattice"
if b_save_figures:
	plt.savefig(s_path + s_file)

plt.show()
tmp = 2

