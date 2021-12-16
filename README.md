# Introduction

This package contains a solver simulating the time evolution of a noisy quantum system of coupled two-level qubits, modeled by a Lindblad master equation.

The code uses matrix-product-state (MPS) and matrix-product-operator (MPO) data structures, implemented in C++ for maximizing performance using multi-threaded computations, and wrapped by a Python layer with an easy-to-use interface and rich plotting features.

The solver performs integration using fixed time steps, of the Lindblad master equation:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial}{\partial t}\rho = -\frac{i}{\hbar}[H,\rho]%2b\mathcal{D}[\rho]">

where <img src="https://render.githubusercontent.com/render/math?math=\rho" style="vertical-align:bottom"> is the density matrix of a quantum system, and
<img src="https://render.githubusercontent.com/render/math?math=[\cdot,\cdot]" align=middle> is the commutator of two operators.

The supported Hamiltonian and dissipator terms of the master equation have time-independent coefficients parametrized as follows;

<img src="https://render.githubusercontent.com/render/math?math={H}/{\hbar} = \sum_{i}\frac{1}{2}\left(h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right) %2b \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right)" align=middle>

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right) %2b \sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right) %2b \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)" style="vertical-align:bottom">

with all the parameters of the model easily specified, and various initial states and observables supported.

For detailed background, information, and examples, the following sections of the documentation can be consulted:

## Table of Contents

* [Installation guide](INSTALL.md)
* [Background and litterature on the method](docs/background.md)
* [Simulator model](docs/MODEL.md)
* [Python interface](docs/API_DOCS.md)
* [C++ Code structure](docs/CODE_STRUCTURE.md)
* [Tutorial: Solving the Lindblad dynamics of a qubit chain](examples/qubit_chain_tutorial.ipynb)
