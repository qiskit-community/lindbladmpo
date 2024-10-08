
# Background and literature on the method

## Basics of the Physics of the problem

This code simulates the time evolution of a quantum system of interacting two-level systems (qubits).

The dynamics is determined by: (1) a Hamiltonian, which corresponds to the unitary part of the time evolution, and (2) dissipative terms, which account for the fact that the system is coupled to an environment.
These two types of terms enter the so-called *Lindblad equation* that determines the evolution of the *density matrix* of the system.
Since we have $N$ qubits the Hilbert space has dimension $2^N$ and the density matrix is $2^N$ by $2^N$ in size. In practice this is a huge dimension unless $N$ is very small, and it therefore prevents a direct brute-force numerical solution of the Lindblad equation.

The present code offers an approximate solution of the problem that can be very accurate for large systems (typically up to $N\sim 100$ or more) if the geometry of the couplings between the qubits is one-dimensional. This approach can also be more efficient than a brute force approach in other geometries.

Two publications in which some simulation results obtained with this code have been described:

https://doi.org/10.1103/PhysRevLett.124.043601

https://doi.org/10.1103/PhysRevB.102.064301

## Matrix product states and matrix-product operators

An MPS is a particular way to encode
a many-body wave-function using a set of matrices. Consider a system made of
$N$ qubits, in a pure state
  ```math
\left|\psi\right\rangle=\sum_{s_1,s_2,\cdots,s_N}\psi\left(s_1,s_2,\cdots,s_N\right)\left|s_1\right\rangle\left|s_2\right\rangle\cdots\left|s_N\right\rangle.
``` 
In this expression the sum runs over the $2^N$ basis states ($s_i\in\\{0,1\\}$) and the wave-function is encoded into the function $\psi:\\{s_i\\}\to\psi\left(s_1,s_2,\cdots,s_N\right).$ An MPS is a state where the wave function is written
  ```math
\psi\left(s_1,s_2,\cdots,s_N\right)={\rm Tr}\left[A^{(s_1)}_1A^{(s_2)}_2\cdots A^{(s_N)}_N\right].
```
where, for each qubit $i$ we have introduced two matrices $A^{(0)}\_i$ and $A^{(1)}\_i$ (for a local Hilbert space of dimension $D$ one needs $D$ matrices $A^{(0)}\_i\cdots A^{(D)}\_i$ for each qubit). These matrices are in general rectangular ($d_i\times d_{i+1}$) and the wave-function is obtained by multiplying them.

What determines the dimensions of the matrices? If the matrices are one-dimensional (scalar) one has  a trivial product state (and all product states can be written this way). On the other hand, if one allows for very large matrices, of size $2^N$, any arbitrary state can be written as an MPS. In fact the MPS representation is really useful when the system has a moderate amount of bipartite entanglement. As a "rule of thumb", to get a good MPS approximation of a given state, each matrix $A^{(s_i)}\_i$, of size $d_{i-1}\times d_{i}$, should have a dimension $d_i$ of the order of $e^{S_{\rm vN}(i)}$, where $S_{\rm vN}(i)$ is the von Neumann entropy of the subsystem $[i+1,\cdots,N]$.


What about mixed states? They can be represented using so-called matrix-product operators (MPO),
  ```math
\rho=\sum_{a_1,a_2,\cdots,a_N}{\rm Tr}\left[M^{(a_1)}_1 M^{(a_2)}_2\cdots M^{(a_N)}_N\right]\sigma^a_1 \otimes \sigma^a_2 \otimes \cdots \sigma^a_N
```
where each $a\_i$ can take four values in $\\{1,x,y,z\\}$, $\sigma^{a}\_i$ is a Pauli matrix or the identity acting on qubit $i$, and we have associated four matrices $M\_i^{(1)}, $M_i^{(x)}$, $M_i^{(y)}$ and $M_i^{(z)}$ to each qubit.


### Bond dimension and entanglement
MPS can allow to store reliably a quantum state with a huge memory gain (compression) when the state is not too entangled. As explained above, 
we should expect that, for a given level of accuracy, the size of the matrices in an MPS will scale as the exponential of the bipartite entanglement entropy (associated to the bi-partition on that 'bond'). In the context of MPS the matrix sizes
$d_i$ are called 'bond dimensions'.

The representation of the pure state in terms of an MPS is therefore all the more efficient as the corresponding pure-state is weakly entangled. By 'more efficient' we mean here that, for a given level of accuracy, the size of the matrices in the MPS will be smaller. Alternatively, if the matrix sizes are fixed, a weakly entangled state will be more precisely represented by an MPS than a highly entangled one.

For a review see for instance: https://doi.org/10.1016/j.aop.2010.09.012

In a parallel way, an operator acting linearly on a many-body state (like a Hamiltonian, an osbservable, or a time-evolution operator) can be encoded as a *matrix-product operator* (MPO). While short-range Hamiltonian in 1D can be expressed exaclty as MPO with a finite bond dimenion which does not depend on the system size, it is not the case in higher dimension or when the MPO encodes a nontricial density matrix.


### DMRG
There are also several powerful alogithms to compute and manipulate quantum states in the form of MPS. The most famous one is the celebrated Density-Matrix Renormalizatio Group (aka DMRG), introduced by S. R. White in 1992 (https://doi.org/10.1103/PhysRevLett.69.2863). It allows to approximate the ground-state of the Hamiltonian of a one-dimensional (1D) system with short-ranged interactions in the form of an MPS. It can also be extended to longer range interactions and to two-dimensional systems, although the system sizes that can be studied in the latter cases are smaller than for 1D systems. The applicability and the efficiency of DMRG again depends on the amount of entanglement present in the targetted states.

We again refer here to the review by U. Schollwöck: https://doi.org/10.1016/j.aop.2010.09.012

## Mixed states

### Vectorization

A mixed state can be viewed as a pure state in some enlarged Hibert space with dimension squared.
This is the so-called *vectorization*, and it is heavily used in this code. For a single qubit a density matrix can be viewed as one vector (= a pure state of some fictitious system) in a space of dimension 4.
In this code a many-body density matrix is considered as a pure state of a (fictitious) system with $(2^N)^2 = 4^N$ states.  In turn, such a pure state is encoded as an MPS (of a system with 4 states per site). The Lindblad super-operator acts linearly on density matrices. Since the present implementation encodes the density matrix as an MPS, the Lindblad super-operator is naturally encoded as an MPO.

This can sometimes be source of confusion: the density matrix is of course an operator acting on the physical Hilbert space of the qubits, but, after vectorization, we interpret it as a pure state, and thus as an MPS.

Some references relevant to the use of MPS and MPO for quantum dissipative systems:

https://doi.org/10.1103/PhysRevLett.93.207204

https://doi.org/10.1103/PhysRevLett.93.207205

http://dx.doi.org/10.1088/1742-5468/2009/02/P02035


### MPS dimension and operator space entanglement entropy

The representation of the many-body density matrix of the system in terms of an MPS is all the more efficient as the corresponding pure state (of the fictitious system) is weakly entangled.  So, it is natural to consider the  von Neumann entanglement entropy of the pure-state (of the fictitious system) obtained from the vectorization of the density matrix (of the real system). We stress that this quantity is not the von Neuman entropy of the real system, it is instead called the Operator Space Entanglement Entropy (OSEE) [https://doi.org/10.1103/PhysRevA.76.032316]. So, for a given target accuracy, the smaller the OSEE the smaller the bond dimension of the MPS (and the fastest the numerical calculations). In practice we instead often fix some maximum bond dimension for the MPS. Then,  the  smaller the OSEE of the physical state the better the MPS approximation will be.

# iTensor Library

The present code is based on the iTensor library, C++ version 3 https://www.itensor.org.
See also the following paper: https://arxiv.org/abs/2007.14822 . The library allows to construct and manipulate MPS and MPO in a simple way and it allows to run DMRG calculations. It also allows to compute the time-evolution of a system where the state is encoded as an MPS and its Hamiltonian is encoded as an MPO.
