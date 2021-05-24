from python.MPOLindbladSolver import MPOLindbladSolver
from qiskit.visualization.gate_map import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Simulation parameters
n_qubits = 5
s_path = "C:/Users/galvz/PycharmProjects/sim_func/"
s_executable = "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator/lindblad.exe"
h_x = np.random.randn(n_qubits).tolist()
h_y = np.random.randn(n_qubits).tolist()
h_z = np.random.randn(n_qubits).tolist()
g_1 = (0.1 * np.random.rand(n_qubits)).tolist()
J = .5
t_final = 10
tau = 1
l_x = 0

# create the parameters dictionary
dict_val = {'tau': tau, 't_final': t_final, 'N': n_qubits,
            'h_x': h_x, 'h_x': h_y, 'h_x': h_z, 'g_1': g_1, 'J': J, 'l_x': l_x,
            'input_file': s_path + "MPO.input",
            'output_file': s_path + "MPO"}
# initialize class - parameters, cygwin_path, simulator_path
sim = MPOLindbladSolver(dict_val, "C:/cygwin64/bin/bash.exe", s_executable)

# execute simulator
sim.solve()

data = np.zeros(shape=(n_qubits, int(t_final/tau) + 1))
i = 0
for key in sim.result['1q']:
    if key[1] == 'Z':
        data[key[0] - 1, i] = sim.result['1q'][key]
        if key[0] == n_qubits:
            i += 1

ax = plt.subplot()
im = ax.imshow(data, interpolation='none')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

ax.set_xlabel('time')
ax.set_ylabel('qubits')
ax.set_title('Pauli Z expectation value as function of time')
plt.show()

