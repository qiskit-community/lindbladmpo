
# Background and litterature on the method

## Basics of the Physics of the problem

This code simulates the time-evolution of a quantum system made of interacting/coupled two-level systems (qubits).

The dynamics is determined by: (1) A Hamiltonian, which corresponds to the unitary part of the time-evolution, and (2) Dissipative terms, which take into account the fact that the system is coupled to an environment.
The two types of terms above enter in a so-called *Lindblad equation* that determines the evolution of the density matrix of the system.
At any given time the state of the system is defined by a (many-body) *density matrix*. Since we have N qubits the Hilbert space has dimension 2^N and the density matrix is 2^N by 2^N in size. In practice this is a huge dimension unless N is not very small, and it therefore prevents a direct brute-force numerical solution of the Lindblad equation.

The present code offers an approximate solution of the problem that can be very accurate for large systems (typically up to N~100 or more) if the geometry of the couplings between the qubits is one-dimensional. This approach can also be more efficient than a brute force approach in other geometries.

Two publications in which some simulation results obtained with this code have been described:

https://doi.org/10.1103/PhysRevLett.124.043601

https://doi.org/10.1103/PhysRevB.102.064301

## Matrix product states and matrix-product operators


### Bond dimension and entanglement
Quantum many-body wavefunctions (i.e. pure states) can be represented as *matrix-product states* (MPS), and this representation often allows to store reliably a quantum state with a huge memory gain (compression). Such a gain is possible when the pure state is not too entangled. A popular measure of the amount of (bipartite) entanglement in a pure state is the von Neumann entropy, also called the entanglement entropy.
Loosely speaking we should expect that, for a given level of accuracy, the size of the matrices in an MPS will scale as the exponential of the bipartite entanglement entropy (associated to the bi-partition on that 'bond'). Note that in the context of MPS the matrix sizes are also called 'bond dimensions'.

The representation of the pure state in terms of an MPS is therefore all the more efficient as the corresponding pure-state is weakly entangled. By 'more efficient' we mean here that, for a given level of accuracy, the size of the matrices in the MPS will be smaller. Alternatively, if the matrix sizes are fixed, a weakly entangled state will be more precisely represented by an MPS than a highly entangled one.

For a review see for instance: https://doi.org/10.1016/j.aop.2010.09.012

In a parallel way, an operator acting linearly on a many-body state (like a Hamiltonian, an osbservable, or a time-evolution operator) can be encoded as a *matrix-product operator* (MPO).

### DMRG
There are also several powerful alogithms to compute and manipulate quantum states in the form of MPS. The most famous one is the celebrated Density-Matrix Renormalizatio Group (aka DMRG), introduced by S. R. White in 1992 (https://doi.org/10.1103/PhysRevLett.69.2863). It allows to approximate the ground-state of the Hamiltonian of a one-dimensional (1D) system with short-ranged interactions in the form of an MPS. It can also be extended to longer range interactions and to two-dimensional systems, although the system sizes that can be studied in the latter cases are smaller than for 1D systems. The applicability and the efficiency of DMRG again depends on the amount of entanglement present in the targetted states.

We again refer here to the review by U. Schollwöck: https://doi.org/10.1016/j.aop.2010.09.012

## MPS, MPO and mixed states

### Vectorization

A mixed state can be viewed as a pure state in some enlarged Hibert space with dimension squared. This is the so-called *vectorization*, and it is heavily used in this code. For a single qubit a density matrix can be viewed as one vector (= a pure state of some fictitious system) in a space of dimension 4. In this code a many-body density matrix is considered as a pure state of a (fictitious) system with (2^N)^2 = 4^N states.  In turn, such a pure state is encoded as an MPS (of a system with 4 states per site). The Lindblad super-operator acts linearly on density matrices. Since the present implementation encodes the density matrix as an MPS, the Lindblad super-operator is naturally encoded as an MPO.

This can sometimes be source of confusion: the density matrix is of course an operator acting on the physical Hilbert space of the qubits, but, after vectorization, we interpret it as a pure state, and thus as an MPS.

Some references relevant to the use of MPS and MPO for quantum dissipative systems:

https://doi.org/10.1103/PhysRevLett.93.207204

https://doi.org/10.1103/PhysRevLett.93.207205

http://dx.doi.org/10.1088/1742-5468/2009/02/P02035


### MPS dimension and operator space entanglement entropy

The representation of the many-body density matrix of the system in terms of an MPS is all the more efficient as the corresponding pure state (of the fictitious system) is weakly entangled.  So, it is natural to consider the  von Neumann entanglement entropy of the pure-state (of the fictitious system) obtained from the vectorization of the density matrix (of the real system). We stress that this quantity is not the von Neuman entropy of the real system, it is instead called the Operator Space Entanglement Entropy (OSEE) [https://doi.org/10.1103/PhysRevA.76.032316]. So, for a given target accuracy, the smaller the OSEE the smaller the bond dimension of the MPS (and the fastest the numerical calculations). In practice we instead often fix some maximum bond dimension for the MPS. Then,  the  smaller the OSEE of the physical state the better the MPS approximation will be.

# iTensor Library

The present code is based on the iTensor library, C++ version 3 https://www.itensor.org.
See also the following papers: https://arxiv.org/abs/2007.14822 . The library allows to construct and manipulate MPS and MPO in a simple way and it allows to run DMRG calculations. It also allows to compute the time-evolution of a system where the state is encoded as an MPS and its Hamiltonian is encoded as an MPO.
