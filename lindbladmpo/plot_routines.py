from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from lindbladmpo.LindbladMPOSolver import LindbladMPOSolver


def prepare_plot_data(solver: LindbladMPOSolver) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
	tau = solver.parameters['tau']
	n_qubits = solver.parameters['N']
	t_final = solver.parameters['t_final']
	# output_step = solver.parameters.get('output_step', 1)
	n_t_steps = int(t_final / (tau * 1)) + 1
	data = np.full(shape = (n_qubits, n_t_steps), dtype = float, fill_value = np.nan)
	i = 0
	for key in solver.result['1q']:
		if key[1] == 'Z':
			data[key[0] - 1, i] = solver.result['1q'][key]
			if key[0] == n_qubits:
				i += 1
	t_steps = np.arange(0, n_t_steps, int(n_t_steps / 10))
	t_ticks = np.round(t_steps * tau, 5)
	qubits = np.asarray(range(n_qubits))
	return data, t_steps, t_ticks, qubits


def plot_space_time(data: np.ndarray, t_steps: np.ndarray, t_ticks: np.ndarray, qubits: np.ndarray,
					ax = None, fontsize = 16, b_save_figures = True, s_file_prefix = ''):
	if ax is None:
		_, ax = plt.subplots(figsize = (14, 9))
	plt.rcParams.update({'font.size': fontsize})
	im = ax.imshow(data, interpolation = 'none', aspect = 'auto')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size = '5%', pad = 0.05)
	plt.colorbar(im, cax = cax)
	# plt.colorbar(im)
	ax.set_xlabel('$t$', fontsize = fontsize)
	ax.set_xticks(t_steps)
	ax.set_xticklabels(t_ticks, fontsize = fontsize)
	ax.set_yticks(qubits)
	ax.set_yticklabels(qubits, fontsize = fontsize)
	ax.set_ylabel('qubits', fontsize = fontsize)
	ax.set_title('$\\langle\sigma_z(t)\\rangle$')
	if b_save_figures:
		plt.savefig(s_file_prefix + '.sigma_z.png')

