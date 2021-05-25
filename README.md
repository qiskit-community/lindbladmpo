
# Background and litterature on the method

## Basics of the Physics of the problem

This code simulates the time-evolution of a quantum system made of interacting/coupled  2-level systems (q-bits)

The dynamics is determined by:
    - Some Hamiltonian, which corresponds to the unitary part of the time-evolution.
    - Some dissipative terms, which take into account the fact that the system is coupled to some environment.
The two types of terms above enter in a so-called *Lindblad equation* which dictates the evolution of the density-matrix of the system.
At any given tiume the state of the system is defined by a (many-body) *density matrix*. Since we have N q-bits the Hilbert space has dimension 2^N and the density matrix is a 2^N * 2^N matrix. In practice this is a huge dimension if N is not very small, and it therefore prevents a direct brute force numerical solution of the Lindblad equation.

The present code offers an approximate solution of the problem which can be  be very accurate for large systems (typically up to N~100 or more) if the geometry of the couplings between the q-bits is one-dimensional. This approach can also be more efficient than a brute force ce approach in other geometries.

Here are two publications in which some simulation results obtained with this code have been described:
https://doi.org/10.1103/PhysRevLett.124.043601
https://doi.org/10.1103/PhysRevB.102.064301

## Matrix product states and matrix-product operators


### Bond dimension and entanglement
Quantum many-body wave-functions (i.e. pure states) can be represented as *matrix-product states* (MPS), and this representation often allows to store reliably a quantum state with a huge memory gain (compression). Such a gain is possible when the pure state is not too entangled. A popular measure of the amount of (bipartite) entanglement in a pure state is the Von Neumann  entropy, also called the entanglement entropy.
Loosely speaking we should expect that, for a given level of accuracy, the size of the matrices in an MPS will scale as the exponential of the bipartite entanglement entropy (associated to the bi-partition on that 'bond'). Note that in the context of MPS the matrix sizes are also called 'bond dimensions'.

The representation of the pure state  in terms of an MPS is therefore all the more efficient as the corresponding pure-state is weakly entangled. By 'more efficient' we mean here that, for a given level of accuracy, the size of the matrices in the MPS will be smaller. Alternatively, if the matrix sizes are fixed, a weakly entangled state will be more precisely represented by an MPS than a highly entangled one.

We again refer here to the review by U. Schollwöck: https://doi.org/10.1016/j.aop.2010.09.012

In a parallel way, an operator acting linearly on a many-body state (like an Hamiltonian, an osbservable, or a time-evolution operator) can be encoded as a *matrix-product operator* (MPO).  

### DMRG
There are also several powerful alogithms to compute and manipulate quantum states in the form of MPS. The most famous one is the celebrated Density-Matrix Renormalizatio Group (aka DMRG), introduced by S. R. White in 1992 (https://doi.org/10.1103/PhysRevLett.69.2863). It allows to approximate the ground-state of the Hamiltonian of a one-dimensional (1D) system with short-ranged interactions in the form of an MPS. It can also be extended to longer range interactions and to two-dimensional systems, although the system sizes that can be studied in the latter cases are smaller than for 1D systems. The applicability and the efficiency of DMRG again depends on the amount of entanglement present in the targetted states.
For a review se for instance https://doi.org/10.1016/j.aop.2010.09.012

## MPS, MPO and mixed states

### Vectorization

A mixed state can be viewed as a pure state in some enlarged Hibert space with dimension squared. This is the so-called *vectorization*, and it is heavily used in this code. For a single q-bit a density matrix can  be viewed as one vector (=pure state of some fictitious system) in a space of dimension 4.  In this code a many-body density matrix is considered as a pure state of a (fictitious) system with (2^N)^2=4^N states.  In turn, such a pure state is encoded as an MPS (of a system with 4 states per site).  The Lindblad super-operator acts linearly on density matrices. Since the present implementation encodes the density matrix as an MPS, the Lindblad super-operator is naturally encoded as an MPO.

This can sometimes be source of confusion: the density matrix is of course an operator acting on the physical Hilbert space of the q-bits, but, after vectorization, we interpret it as a pure state, and thus as an MPS. 

Some references relevant to the use of MPS and MPO for quantum dissipative systems: https://doi.org/10.1103/PhysRevLett.93.207204, https://doi.org/10.1103/PhysRevLett.93.207205 and http://dx.doi.org/10.1088/1742-5468/2009/02/P02035


### MPS dimension and operator space entanglement entropy

The representation of the many-body density matrix of the system in terms of an MPS is all the more efficient as the corresponding pure-state (of the fictitious system) is weakly entangled.  So, it is natural to consider the  Von Neumann entanglement entropy of the pure-state (of the fictitious system) obtained from the vectorization of the density-matrix (of the real system). We stress that this quantity is not the Von Neuman entropy of the real system, it is instead called the Operator Space Entanglement Entropy (OSEE) [https://doi.org/10.1103/PhysRevA.76.032316]. So, for a given target accuracy, the smaller the OSEE the smaller the bond dimension of the MPS (and the fastest the numerical calculations). In practice we instead often fix some maximum bond dimension for the MPS. Then,  the  smaller the OSEE of the physical state the better the MPS approximation will be.

# iTensor Library

The present code is based on the iTensor library, version 3 https://www.itensor.org.
See also the following papers: https://arxiv.org/abs/2007.14822 . The library allows to construct and manipulate MPS and MPO in a simple way and it allows to run DMRG calculations. It also allows to compute the time-evolution of a system where the state is encoded as an MPS and its Hamiltonian is encoded as an MPO.

# A few words about the code structure

* The files `Pauli.cc` and `Pauli.h` contain the implementation of the vectorization for spin-1/2 systems, that is the representation of the many-body density matrix using an MPS. It is specific to 2-level systems (q-bits), but not specific the special form of the Lindbadian of the specific model to be studied. It makes use of a few low-level routines of the iTensor library. Remark: The identity matrix plus the three Pauli matrices form a convenient basis of the the space of density matricxes of one q-bit, hence the name of these two files. 

* `lindblac.cc`: contains the main() function. This is where the initial state is constructed (the way the initial state is defined depends on the input parameters), and where the loop over the time steps is defined. It is also where the observable of interests are computed, and the input and outputs (simple prints to the standard output, or to files) are handled.

Possible initial states: 1) A pure (product) state where all the spins point in a weel-defined direction 2) An inifinite-temperature (mixed) state, where the density matrix is proportionnal to the identity matrix 3) A mixed state which is read from the disk. Such a state may have been written to disk at the end of some previous simulation.

* `Mylindbladian.h`: This is the place where the Lindbladian super-operator of the specific model to be simulated is defined. It takes the form of an `autoMPO` object of the iTensor library. Thanks to the use of the iTensor library, the terms in the Lindbladian can be defined in a simple way (using spin-1/2 operators and product of such operators), almost as if one was writing them with pen and paper.

* `ModelParameters.h`: defines the parameters which are specific to the model to be simulated. These parameters are couples of the type "name" / "value" (values are floating point numbers or integers). Their values can be defined in the command-line when calling the executable. Example `./lindblad param1 val1 0.2 val2 10`. Each parameter must have a default value, so that a parameter not specified in the command line will take its default value. The default values are defined in this same file. The same arguments can also be specified in some input file. Example:  `./lindblad inputfile in.txt` where the file `in.txt` contains lines of the form variable1 = value_of_variable_1 etc. 

* `SimulationParameters.h`: defines the parameters which are specific to the MPS/MPO simulation method, and in particular about the wanted level of  accuracy (maximum truncation, maximum bond dimension, etc.)    

* `SimpleSquareLattice.h`: contains a basic class to encode the spatial geometry of the system (a lattice). In the present version it can handle a simple one-dimensional chain, or a square lattice with cylindrical boundary conditions. Similar classes could be created to handle other geometries and other lattices. 

* `TimeEvolution.h` and `TimeEvolution.cc`: Contain the `TimeEvolver` class. Such class stores the parameters associated to a single-time-step evolution of the density matrix. An important parameter is for instance the value/length `tau` of one  time step. A `TimeEvolver` also contains other parameters associated to the approximations (truncations, etc.) to be made when applying such the time evolution operator (which is an MPO) to a given density matrix. The `TimeEvolver` class is independent of the details of the specific model to be studied. The actual time-evolution is coded in `TimeEvolution.cc`: the method `evolve` takes a density matrix as an input [it is an iTensor MPS], an updates it 'in place' by the evolved one. Different Trotter orders are available. At order o=2 the error made at each time state is O(tau^3). At order o=3 the error made at each time state is O(tau^4.). At order o=4 the error made at each time state is O(tau^5).

* `io_util.h`: small and  simple methods for inputs and outputs (not specific to this type of simulations). It contains the class `Parameters`, which is used to handle a set of input parameters (each one being a pair "name" / "value") defined from a acommand line.

* `mps_mpo_utils.h` and `mps_mpo_utils.h`: basic but useful general methods for MPS and MPO (not specific to density matrices nor dissipative systems).
