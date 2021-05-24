from python.MPOLindbladSolver import MPOLindbladSolver
from qiskit.visualization.gate_map import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Simulation parameters
n_qubits = 17
# s_path = "C:/Users/galvz/PycharmProjects/sim_func/"
# s_executable = "/cygdrive/c/Users/galvz/AppData/Roaming/SPB_Data/Lindbladian-MPO-simulator/lindblad.exe"
s_path = "C:/temp/"
s_executable = "/cygdrive/c/Users/HaggaiLanda/gitprojects/Lindbladian-MPO-simulator/lindblad.exe"
b_save_figures = True
fontsize = 16

h_x = (0. * np.random.randn(n_qubits)).tolist()
h_x[int(n_qubits / 2)] = 2.
h_y = (0. * np.random.randn(n_qubits)).tolist()
h_z = (0. * np.random.randn(n_qubits)).tolist()
g_1 = (0.01 * np.random.rand(n_qubits)).tolist()
J = 4
t_final = 1
tau = .02
l_x = n_qubits

# create the parameters dictionary
dict_val = {'tau': tau, 't_final': t_final, 'N': n_qubits,
            'h_x': h_x, 'h_z': h_z, 'g_1': g_1, 'J': J, 'l_x': l_x,
            'input_file': s_path + "MPO.input",
            'output_file': s_path + "MPO"}
# initialize class - parameters, cygwin_path, simulator_path
sim = MPOLindbladSolver(dict_val, "C:/cygwin64/bin/bash.exe", s_executable)

# execute simulator
sim.solve()

n_t_steps = int(t_final/tau) + 1
data = np.zeros(shape=(n_qubits, n_t_steps))
i = 0
for key in sim.result['1q']:
    if key[1] == 'Z':
        data[key[0] - 1, i] = sim.result['1q'][key]
        if key[0] == n_qubits:
            i += 1
t_steps = np.arange(0, n_t_steps, int(n_t_steps / 10))
t_ticks = np.round(t_steps * tau, 5)
r_qubits = range(n_qubits)

fig, ax = plt.subplots(figsize = (14, 8))
plt.rcParams.update({'font.size': fontsize})
im = ax.imshow(data, interpolation='none', aspect = 'auto')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
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

