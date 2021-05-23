from python.MPOLindbladSolver import MPOLindbladSolver
from qiskit.visualization.gate_map import *

n_qubits = 5
s_path = "C:/temp/"
s_executable = "/cygdrive/c/Users/HaggaiLanda/gitprojects/Lindbladian-MPO-simulator/lindblad.exe"
h_x = np.random.randn(n_qubits).tolist()
h_y = np.random.randn(n_qubits).tolist()
h_z = np.random.randn(n_qubits).tolist()
g_1 = (0.1 * np.random.rand(n_qubits)).tolist()
J = .5
t_final = 10
tau = 0.02

# create the parameters dictionary
dict_val = {'tau': tau, 't_final': t_final, 'N': n_qubits,
            'h_x': h_x, 'h_x': h_y, 'h_x': h_z, 'g_1': g_1, 'J': J,
            'input_file': s_path + "MPO.input",
            'output_file': s_path + "MPO"}
# initialize class - parameters, cygwin_path, simulator_path
sim = MPOLindbladSolver(dict_val, "C:/cygwin64/bin/bash.exe", s_executable)

# execute simulator
sim.solve()
