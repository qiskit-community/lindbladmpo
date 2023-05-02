# Introduction

**This package contains a solver simulating the time evolution of an open quantum system of coupled two-level qubits, modeled by a Lindblad master equation.**

The code uses matrix-product-state (MPS) and matrix-product-operator (MPO) data structures, implemented in C++ (using the [ITensor](https://itensor.org/) library) for maximizing the performance with multithreaded computations, and wrapped by a Python layer with an easy-to-use interface and rich plotting features. The solver integrates in fixed time steps the Markovian master equation for the density matrix,

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial}{\partial t}\rho = -\frac{i}{\hbar}[H,\rho]%2b\mathcal{D}[\rho]">

The solver supports various initial states and observables, with the Hamiltonian and dissipator terms having time-independent coefficients;

<img src="https://render.githubusercontent.com/render/math?math={H}/{\hbar} = \sum_{i}\frac{1}{2}\left(h_{z,i}\sigma_i^z  %2b h_{x,i}\sigma_i^x %2b h_{y,i}\sigma_i^y\right) %2b \frac{1}{2}\sum_{ i}^N\sum_{ j\neq i}^N \left(J_{ij}\sigma^x_i \sigma^x_{j} %2b J_{ij}\sigma^y_i \sigma^y_{j} %2b J_{ij}^z \sigma^z_i \sigma^z_{j}\right)," align=middle>

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{D}[\rho] = \sum_i g_{0,i}\left(\sigma_i^%2b \rho\sigma_i^- - \frac{1}{2} \{\sigma_i^- \sigma_i^%2b,\rho\}\right) %2b \sum_i g_{1,i}\left( \sigma_i^-\rho \sigma_i^{%2b}-\frac{1}{2}\left\{\sigma_i^{%2b}\sigma_i^-,\rho\right\}\right) %2b \sum_i g_{2,i} \left(\sigma_i^z \rho\sigma_i^z - \rho\right)." style="vertical-align:bottom">

**How to cite:** The following paper introduces the solver and presents a research performed using this solver:

H. Landa and G. Misguich, *Nonlocal correlations in noisy multiqubit systems simulated using matrix product operators*,
[SciPost Phys. Core **6**, 037 (2023)](https://scipost.org/SciPostPhysCore.6.2.037).

**Verifying numerical simulations:** The above paper also explains in detail how to construct and verify numerical results with multiqubit simulations. We recommend users of the package to familiarize themselves with different approaches to verify their numerical simulations, such as those explained in this paper. The sources used to perform this research are included in the detailed [examples](docs/examples.md).

Note that the C++ binaries of the solver must be **built locally** - see the [installation guide](INSTALL.md).

For more information, tutorials, and detailed examples check the documentation:

## Table of Contents

* [Installation guide](INSTALL.md)
* [Tutorial: Solving the Lindblad dynamics of a qubit ring](docs/qubit_chain_tutorial.ipynb) (go here for a **basic example**)
* [Background and literature on the solver method](docs/background.md)
* [Supported model](docs/dynamical_model.md)
* [Python interface](docs/API_DOCS.md) (**API documentation** of the Python input, output, and call interface)
* [Simulation examples](docs/examples.md)  (**advanced** Python examples for research using the solver)
* [C++ solver interface](docs/cpp_solver_interface.md) (only relevant for an advanced, low-level usage)
* [C++ code structure](docs/cpp_code_structure.md) (only relevant for contributing to the C++ code)
