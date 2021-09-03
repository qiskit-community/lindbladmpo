# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
import matplotlib.pyplot as plt

from lindbladmpo.LindbladMPOSolver import *
from lindbladmpo.plot_routines import *
import os.path
import csv
import pandas as pd

s_output_path = os.path.abspath('./output') + '/'
if not os.path.exists(s_output_path):
    os.mkdir(s_output_path)
s_file_prefix = "chain"

# Simulation parameters
rand_seed = 1
n_qubits = 7
b_unique_id = True
b_save_to_db = True
b_save_figures = False
fontsize = 22
h_z_half_width = 0.
oneq_components = ['Z', 'Y']
twoq_components = ['ZZ', 'XX', 'XY']

np.random.seed(rand_seed)
h_x = np.zeros(n_qubits, float)
h_x[int(n_qubits / 2)] = 5.

h_z = 0
g_0 = 0
J = 1
t_final = 1
tau = 0.1

db_default_params = {'id': 0, 'mapping': 'A', 'topology': 'chain', 'solver': 'MPO', 'N': 7, 'tau': 1,
                     't_final': 1, 'cut_off_rho': 1e-9, 'max_dim_rho': 500, 'h_x': 0, 'h_z': 0, 'g_2': 0}

# Create the parameters dictionary
solver_params = {'tau': tau, 't_final': t_final, 'cut_off_rho': 1e-16,
                 'max_dim_rho': 500, 'N': n_qubits, 'b_unique_id': b_unique_id,
                 'h_x': h_x, 'h_z': h_z, 'g_0': g_0, 'J': J, 'l_x': n_qubits, 'l_y': 1,
                 'output_files_prefix': s_output_path + s_file_prefix, '1q_components': oneq_components,
                 '2q_components': twoq_components}
# Initialize class arguments - the parameters, cygwin path, and MPO executable path
s_cygwin_path = None
s_solver_path = None
solver = LindbladMPOSolver(solver_params, s_cygwin_path, s_solver_path)
# Execute simulator
#solver.solve()

header = db_default_params.keys()
dict_to_db = {}
print(os.path.isfile("database.csv"))
if not os.path.isfile("database.csv"):
    with open('database.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        f.close()

# open the file in the write mode
if b_save_to_db:
    df = pd.read_csv("database.csv")
    for key in db_default_params.keys():
        if key in solver_params.keys():
            dict_to_db[key] = [solver_params[key]]
        else:
            dict_to_db[key] = [db_default_params[key]]
    print(dict_to_db)
    data_to_db = pd.DataFrame(dict_to_db)
    print(data_to_db)
    df = df.append(data_to_db)
    df.to_csv('database.csv', index=False)

# Prepare plot data
for obs in oneq_components:
    data, t_steps, t_ticks, qubits = prepare_data_chain(solver, obs)
    for qubit in range(np.shape(data)[0]):
        plt.plot(t_steps, data[qubit])
    plt.legend(qubits)
    plt.title('<' + obs + '>')
    plt.show()

'''
print(data)
print(t_steps)
print(t_ticks)
print(qubits)
'''

# And plot one figure. By default it will be saved to a file.
# plot_space_time(data, t_steps, t_ticks, qubits, fontsize = fontsize, s_file_prefix = s_output_path + s_file_prefix,
# b_save_figures = b_save_figures)
