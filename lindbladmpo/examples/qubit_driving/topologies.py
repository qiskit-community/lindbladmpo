# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Definitions of the topologies for the example research projects.
"""
from typing import Tuple
import matplotlib.pyplot as plt
import numpy as np


coupling_maps = {}
"""A dictionary with coupling maps (lists of bonds) for different qubit configurations and topologies."""

qubit_coordinates = {}
"""A dictionary with coordinates lists positioning qubits in the plane for plotting configurations,
corresponding to those defined in the `coupling_maps` dictionary."""

h_z_patterns = {}
"""A dictionary with a pattern of the relative strength of h_z terms of a qubit configuration."""


def _create_ring_A(n_qubits: int, i_offset: int = 0) -> Tuple[list, list, list]:
	"""Generate a ring topology of type A - indicating a ladder-like ordering of the qubits.

		Args:
			n_qubits: The number of qubits composing the ring.
			i_offset: The index offset of the first qubit.
		Returns:
			A tuple with the coupling map, qubit coordinates in the plane (with the first qubit placed
				at [1, 1]), and a 0-1 pattern array with alternating values for neighboring qubits.
	"""
	c_map = []
	q_coordinates = []
	h_z_pat = []
	c_map.extend([[i_offset, i_offset + 1], [i_offset, i_offset + 2]])
	for i in range(i_offset + 1, i_offset + n_qubits - 2):
		c_map.append([i, i + 2])
	c_map.append([i_offset + n_qubits - 2, i_offset + n_qubits - 1])
	h_z_pat.append(0)
	q_coordinates.append([1, 1])
	h_z_ = 1
	for i in range(1, n_qubits - 2, 2):
		h_z_pat.extend([h_z_, h_z_])
		q_coordinates.extend([[2 + int(i / 2), 0], [2 + int(i / 2), 2]])
		h_z_ = 0 if h_z_ == 1 else 1
	h_z_pat.append(h_z_)
	q_coordinates.append([1 + int(n_qubits / 2), 1])
	return c_map, q_coordinates, h_z_pat


# Add a 1D chain topology for an odd number of qubits (3 to 29 qubits).
# The chain topology entries have keys of the form 'N.chain.M' where N is the number of qubits,
# and M indicates that the middle qubit has h_z amplitude 0.
for n_qubits in range(3, 33, 2):
	c_map = []
	q_coordinates = []
	h_z_pat = []
	for i in range(n_qubits - 1):
		c_map.append([i, i + 1])
	driven_qubit_is_odd = 1 if (n_qubits - 1) % 4 != 0 else 0
	for i in range(n_qubits):
		q_coordinates.append([0, i])
		h_z_pat.append(0 if i % 2 == driven_qubit_is_odd else 1)
	s_key = f'{n_qubits}.chain.M'
	coupling_maps[s_key] = c_map
	qubit_coordinates[s_key] = q_coordinates
	h_z_patterns[s_key] = h_z_pat


# Add a 1D chain topology for all numbers of qubits (3 to 61 qubits).
# The chain topology entries have keys of the form 'N.chain.E' where N is the number of qubits,
# and E indicates that the left edge qubit (indexed 0) has h_z amplitude 0.
for n_qubits in range(2, 63):
	c_map = []
	q_coordinates = []
	h_z_pat = []
	for i in range(n_qubits - 1):
		c_map.append([i, i + 1])
	for i in range(n_qubits):
		q_coordinates.append([0, i])
		h_z_pat.append(0 if i % 2 == 0 else 1)
	s_key = f'{n_qubits}.chain.E'
	coupling_maps[s_key] = c_map
	qubit_coordinates[s_key] = q_coordinates
	h_z_patterns[s_key] = h_z_pat


# We add ring topologies for 4 to 32 qubits, with keys in the form 'N.ring.A' where N is the number
# of qubits, and A is the mpo_ordering - indicating a ladder-like ordering of the qubits.
for n_qubits in range(4, 34, 2):
	c_map, q_coordinates, h_z_pat = _create_ring_A(n_qubits)
	s_key = f'{n_qubits}.ring.A'
	coupling_maps[s_key] = c_map
	qubit_coordinates[s_key] = q_coordinates
	h_z_patterns[s_key] = h_z_pat


# We add ring topologies for 4 to 32 qubits, with keys in the form 'N.ring.B' where N is the number of
# qubits, and B is the mpo_ordering - indicating qubits ordered to have one large jump at the last one.
for n_qubits in range(4, 34, 2):
	c_map = []
	q_coordinates = []
	h_z_pat = []
	for i in range(n_qubits - 1):
		c_map.append([i, i + 1])
	c_map.append([n_qubits - 1, 0])
	for i in range(n_qubits):
		h_z_pat.append(0 if i % 2 == 0 else 1)
	q_coordinates.append([1, 1])
	N_by_2 = int(n_qubits / 2)
	for i in range(2, N_by_2 + 1):
		q_coordinates.append([i, 0])
	q_coordinates.append([N_by_2 + 1, 1])
	for i in range(N_by_2, 1, -1):
		q_coordinates.append([i, 2])
	s_key = f'{n_qubits}.ring.B'
	coupling_maps[s_key] = c_map
	qubit_coordinates[s_key] = q_coordinates
	h_z_patterns[s_key] = h_z_pat


# We add plaquette topologies for 6 to 42 qubits, with keys in the form 'N.plaquette.A' where N is
# the number of qubits, and A is the mpo_ordering - indicating a ladder-like ordering of the qubits.
for n_qubits in range(6, 44, 2):
	c_map1, q_coordinates1, h_z_pat1 = _create_ring_A(n_qubits - 2, 1)
	c_map = [[0, 1]]
	c_map.extend(c_map1)
	c_map.append([n_qubits - 2, n_qubits - 1])
	q_coordinates = [[0, 1]]
	q_coordinates.extend(q_coordinates1)
	q_coord_end = q_coordinates1[-1]
	q_coordinates.append([q_coord_end[0] + 1, q_coord_end[1]])
	h_z_pat = [0]
	for h_z in h_z_pat1:
		h_z_pat.append(1 if h_z == 0 else 0)
	h_z_pat.append(0)
	s_key = f'{n_qubits}.plaquette.A'
	coupling_maps[s_key] = c_map
	qubit_coordinates[s_key] = q_coordinates
	h_z_patterns[s_key] = h_z_pat


# We add some plaquette topologies, with keys in the form 'N.plaquette.B' where N is the number of
# qubits, and B is the mpo_ordering - indicating qubits ordered to have one large jump at the last one.

coupling_maps['10.plaquette.B'] = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
								   [5, 6], [5, 7], [7, 8], [8, 9], [9, 1]]
qubit_coordinates['10.plaquette.B'] = [[0, 1], [1, 1], [2, 0], [3, 0],
									   [4, 0], [5, 1], [6, 1], [4, 2],
									   [3, 2], [2, 2]]
h_z_patterns['10.plaquette.B'] = [0, 1, 0, 1, 0, 1, 0, 0, 1, 0]

coupling_maps['12.plaquette.B'] = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5],
								   [5, 6], [6, 7], [6, 8], [8, 9], [9, 10],
								   [10, 11], [11, 1]]
qubit_coordinates['12.plaquette.B'] = [[0, 1], [1, 1], [2, 0], [3, 0],
									   [4, 0], [5, 0], [6, 1], [7, 1],
									   [5, 2], [4, 2], [3, 2], [2, 2]]
h_z_patterns['12.plaquette.B'] = [0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0]


# Below is the connectivity map and outline of the IBM Quantum Falcon devices

coupling_maps['27.falcon'] = [[0, 1], [1, 0], [1, 2], [1, 4], [2, 1], [2, 3], [3, 2], [3, 5],
							  [4, 1], [4, 7], [5, 3], [5, 8], [6, 7], [7, 4], [7, 6], [7, 10],
							  [8, 5], [8, 9], [8, 11], [9, 8], [10, 7], [10, 12], [11, 8], [11, 14],
							  [12, 10], [12, 13], [12, 15], [13, 12], [13, 14], [14, 11], [14, 13],
							  [14, 16], [15, 12], [15, 18], [16, 14], [16, 19], [17, 18], [18, 15],
							  [18, 17], [18, 21], [19, 16], [19, 20], [19, 22], [20, 19], [21, 18],
							  [21, 23], [22, 19], [22, 25], [23, 21], [23, 24], [24, 23], [24, 25],
							  [25, 22], [25, 24], [25, 26], [26, 25]]
qubit_coordinates['27.falcon'] = [[1, 0], [1, 1], [2, 1], [3, 1], [1, 2], [3, 2], [0, 3], [1, 3],
	[3, 3], [4, 3], [1, 4], [3, 4], [1, 5], [2, 5], [3, 5], [1, 6], [3, 6], [0, 7], [1, 7],
	[3, 7], [4, 7], [1, 8], [3, 8], [1, 9], [2, 9], [3, 9], [3, 10]]
h_z_patterns['27.falcon'] = [0] * 27


def plot_topology(N: int, topology: str, s_coupling_map: str, b_save_figures: bool,
				  s_file_prefix: str, b_transpose_plot = False, b_alternating_qubits = False):
	qubit_color = ["#648fff"] * N
	h_z_pattern = h_z_patterns[s_coupling_map]
	coupling_map = coupling_maps[s_coupling_map]
	qubit_coord = qubit_coordinates[s_coupling_map]
	if b_alternating_qubits:
		s_alternating = '.alternating'
		for i in np.nonzero(h_z_pattern)[0]:
			qubit_color[i] = "#ff6f64"
	else:
		s_alternating = ''

	if topology == 'plaquette' or topology == 'ring':
		figsize = (4, 7)
	else:
		figsize = (8, 2)
	q_coord = []
	if b_transpose_plot:
		figsize = (figsize[1], figsize[0])
		for ll in qubit_coord:
			q_coord.append([ll[1], ll[0]])
	else:
		q_coord = qubit_coord

	try:
		from qiskit.visualization.gate_map import plot_coupling_map
		fig = plot_coupling_map(num_qubits = N, qubit_coordinates = q_coord,
							coupling_map = coupling_map, figsize = figsize, qubit_color = qubit_color)
		if b_save_figures:
			plt.savefig(s_file_prefix + s_alternating + ".png")
		plt.draw()
		plt.pause(0.1)
		plt.show(block = False)
	except Exception as e:
		print(str(e))
