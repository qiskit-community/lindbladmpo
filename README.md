# Introduction

This module simulates the time-evolution of a quantum system made of interacting/coupled two-level systems (qubits).
The code uses matrix-product states (MPS) and matrix-product operator (MPO) structures to efficiently simulate an approximate solution to the evolution of a system containing many qubits.

The simulated system Hamiltonian and Dissipator are in the formats of:

Hamiltonian:
<img src="https://render.githubusercontent.com/render/math?math=\frac{H}{\hbar} = \sum_{i}\frac{1}{2}\left[h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right] %2b \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right)">

Dissipator:
<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right) %2b \sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right) %2b \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)">

while all the coefficients of the model are controled by the user.

for more information see the links below.

# Table of Contents

* [Installation guide](docs/INSTALL.md)
* [Background and litterature on the method](docs/background.md)
* [Simulator's model](docs/MODEL.md)
* [Python interface](docs/API_DOCS.md)
* [Code structure](docs/CODE_STRUCTURE.md)
* [Tutorial: Solving the Lindblad dynamics of a qubit chain] (examples/qubit_chain_tutorial.ipynb)
