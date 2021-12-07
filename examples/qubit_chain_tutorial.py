# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np
from lindbladmpo.LindbladMPOSolver import *
from lindbladmpo.plot_routines import *

s_output_path = os.path.abspath('./output') + '/'
if not os.path.exists(s_output_path):
	os.mkdir(s_output_path)
s_file_prefix = "chain_tutorial3"

####################
# model parameters #
####################

# basic parameters
N = 6
t_init = 0 # this is also the default value
t_final = 2
tau = 0.01
output_files_prefix = s_output_path + s_file_prefix


nu_z = 4
nu_x = 1
XY_interaction = 4
# Hamiltonian coefficients (remember that all the Hamiltonian is multiplied by 0.5)
h_z = 2 * np.pi * nu_z
h_x = 2 * np.pi * nu_x
J = 2 * np.pi * XY_interaction
J_z = 0 # this is also the default value

# Dissipation coefficients
g_0 = 4

# Lattice
# We want to model a ring consisting N qubits
l_x = N
l_y = 1 # this is also the default value
b_periodic_x = True

# Observables
one_qubit_observables = ['x', 'y', 'z']
two_qubits_observables = []

# Create the parameters dictionary
solver_params = {'N': N, 't_init': t_init, 't_final': t_final,  'tau': tau, 'output_files_prefix': output_files_prefix,
				 'h_x': h_x, 'h_z': h_z, 'J': J, 'J_z': J_z, 'g_0': g_0,
                 'l_x': l_x, 'l_y': l_y, 'b_periodic_x': b_periodic_x,
				 '1q_components': one_qubit_observables, '2q_components': two_qubits_observables}

# Initialize the class
solver = LindbladMPOSolver(solver_params)
# Execute simulator
solver.solve()
# solver.result = solver.load_output(solver.s_output_path)

x_expectation_values = np.array([[solver.result['obs-1q'][('x', (i,))][1][t] for t in range(len(solver.result['obs-1q'][('x', (i,))][0]))] for i in range(N)])
y_expectation_values = np.array([[solver.result['obs-1q'][('y', (i,))][1][t] for t in range(len(solver.result['obs-1q'][('y', (i,))][0]))] for i in range(N)])
z_expectation_values = np.array([[solver.result['obs-1q'][('z', (i,))][1][t] for t in range(len(solver.result['obs-1q'][('z', (i,))][0]))] for i in range(N)])

# the expectation values are tuples in the format of ([times], [values])
# sums = []
# for i in range(N):
#     cur_sum = 0
#     for j in range(len(x_expectation_values[i][0])):
#         cur_sum += x_expectation_values[i][1][j]
#     sums.append(cur_sum)

n_steps = x_expectation_values.shape[1]
x_mean = [np.mean(x_expectation_values[:, i]) for i in range(n_steps)]
y_mean = [np.mean(y_expectation_values[:, i]) for i in range(n_steps)]
z_mean = [np.mean(z_expectation_values[:, i]) for i in range(n_steps)]
print(x_mean)

from qiskit.visualization import plot_bloch_vector
import matplotlib.pyplot as plt

fontsize = 16

t_eval = np.linspace(t_init, t_final, n_steps)

_, ax = plt.subplots(figsize = (10, 6))
plt.rcParams.update({'font.size': fontsize})
plt.plot(t_eval, x_mean, label = '$ N^{-1}\sum_i \\langle X_i \\rangle$')
plt.plot(t_eval, y_mean, label = '$ N^{-1}\sum_i \\langle Y_i \\rangle$')
plt.plot(t_eval, z_mean, label = '$ N^{-1}\sum_i \\langle Z_i \\rangle$')
plt.legend(fontsize = fontsize)
ax.set_xlabel('$t$', fontsize = fontsize)
ax.set_title('Mean Bloch vector vs. $t$', fontsize = fontsize)

plot_bloch_vector([x_mean[-1], y_mean[-1], z_mean[-1]], 'Mean Bloch vector at $t = {t_eval[-1]}$')

print(y_mean)