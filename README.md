# Introduction

This package contains a solver simulating the time-evolution of a noisy quantum system of coupled two-level qubits, modeled by a Lindblad master equation.
The code uses matrix-product state (MPS) and matrix-product operator (MPO) structures to efficiently simulate an approximate solution to the evolution of a system containing many qubits.

The solver is implemented in C++ for maximizing performance using multi threaded calculations, and wrapped by a Python layer with rich plotting features and an easy-to-use interface.

The supported Hamiltonian and the Lindbladian dissipators of the master equation have time-independent coefficients parametrized as follows;

<img src="https://render.githubusercontent.com/render/math?math={H}/{\hbar} = \sum_{i}\frac{1}{2}\left(h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right) %2b \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right)" align=middle>

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right) %2b \sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right) %2b \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)" style="vertical-align:bottom">

with all the parameters of the model easily specified, and various initial states and observables supported.

For detailed information browse the documentation using the links below.

# Table of Contents

* [Installation guide](docs/INSTALL.md)
* [Background and litterature on the method](docs/background.md)
* [Simulator model](docs/MODEL.md)
* [Python interface](docs/API_DOCS.md)
* [Code structure](docs/CODE_STRUCTURE.md)
* [Tutorial: Solving the Lindblad dynamics of a qubit chain](examples/qubit_chain_tutorial.ipynb)
