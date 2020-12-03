
# Background and litterature on the method

## Basics of the Physics of the problem

This code simulates the time-evolution of a quantum system made of interacting/coupled  2-level systems (q-bits)

The dynamics is determined by:
    - Some Hamiltonian, which corresponds to the unitary part of the time-evolution.
    - Some dissipative terms, which take into account the fact that the system is coupled to some environment.
The two types of terms above enter in a so-called *Lindblad equation* which dicates the evolution of the density-matrix of the system.
At any given tiume the state of the system is defined by a (many-body) *density matrix*. Since we have N q-bits the Hilbert space has dimension 2^N and the density matrix is a 2^N * 2^N matrix. In practice this is a huge dimension if N is not very small, and it therefore prevents a direct brute force numerical solution of the Lindblad equation.

The present code offers an approximate solution of the problem which can be  be very accurate for large systems (typically up to N~100 or more) if the geometry of the couplings between the q-bits is one-dimensional. This approach can also be more efficient than a brute force ce approach in other geometries.

Here are two publications in which some simulation results obtained with this code have been described:
https://doi.org/10.1103/PhysRevLett.124.043601
https://doi.org/10.1103/PhysRevB.102.064301

## Matrix product states and matrix-product operators

Quantum many-body wave-functions (i.e. pure states) can be represented as *matrix-product states* (MPS), and this representation often allows to store reliably a quantum state with a huge memory gain (compression). There are also several powerful alogithms to compute and manipulate quantum states in the form of MPS. The most famous one is probably the celebrated Density-Matrix Renormalizatio Group (aka DMRG https://doi.org/10.1103/PhysRevLett.69.2863, https://doi.org/10.1016/j.aop.2010.09.012) In the same way an operator acting linearly on a many-body state can be encoded as a *matrix-product operator* (MPO).  

Some references relevant to the use of MPS and MPO for quantum dissipative systems: https://doi.org/10.1103/PhysRevLett.93.207204, https://doi.org/10.1103/PhysRevLett.93.207205 and http://dx.doi.org/10.1088/1742-5468/2009/02/P02035

## Vectorization

A mixed state can be viewed as a pure state in some enlarged Hibert space with dimension squared. This is the so-called *vectorization*, and it is heavily used in this code. For a single q-bit a density matrix can  be viewed as one vector (=pure state) in a space of dimension 4.  In this code a many-body density matrix is considered as a pure state (or a wave-function) in some enlarged Hilbert space of a system with (2^N)^2=4^N states per site.  In turn, such a pure state is encoded as an MPS (of a system with 4 states per site).  The Lindblad super-operator acts linearly on density matrices. Since the present implementation
encodes the density matrix as an MPS, the Lindblad super-operator is naturally encoded as an MPO.
This can sometimes be source of confusion: the density matrix is of course an operator acting on the physical Hilbert space of the q-bits, but, after vectorization, we interpret it as a pure state and thus an MPS. 

## iTensor Library

The present code is based on the iTensor library, version 3 https://www.itensor.org.
See also the following papers: https://arxiv.org/abs/2007.14822

# A few words about the code structure

* The files `Pauli.cc` and `Pauli.h` contain the core of the method, i.e. the representation of the many-body density matrix using an MPS (vectorization). It is specific to 2-level systems (q-bits), but not specific the special form of the Lindbadian of the specific model to be studied. It makes use of a few low-level routines of the iTensor library. Remark: The identity matrix plus the three Pauli matrices form a convenient basis of the the space of density matricxes of one q-bit, hence the name of these two files. 

* `lindblac.cc`: contains the main() function. This is where the initial state is constructed (from the input parameters), and where the loop over the time steps is defined. It is also where the observable of interests are computed, and the input and outputs are handled 

* `Mylindbladian.h`: This is the place where the Lindbladian super operator of the specific model to be simulated is defined. It takes the form of an `autoMPO` object of the iTensor library.

* `ModelParameters.h`: defines the parameters which are specific to the model to be simulated. These parameters are couples of the type "name" / "value" (values are floating point numbers or integers). Their values can be defined in the command-line when calling the executable. Example `./lindblad param1 val1 0.2 val2 10`. Each parameter must have a default value, so that a parameter not specified in the command line will take its default value. The default values are defined in this same file.

* `SimulationParameters.h`: defines the parameters which are specific to the MPS/MPOsimulation method, and in particular about the level of wanted accuracy (maximum truncation, maximum bond dimension, etc.)

* `SimpleSquareLattice.h`: contains a basic class to encode the spatial geometry of the system. In the present version it can handle a simple chain, or a square lattice with cylindrical boundary conditions. Similar classes could be created to handle other geometries.

* `TimeEvolution.h` and `TimeEvolution.cc`: Contain the `TimeEvolver` class. Such class stores the parameters associated to a 1-time-step evolution of the density matrix. The main parameter is for instance the value/length `tau` of one  time step. It also contains other parameters associated to the approximations (truncations, etc.) to be made when applying such a time evolution to a given density matrix. This implementation of this class is independent of the details of the specific model to be studied. The actual time-evolution is coded in `TimeEvolution.cc`: the method `evolve` takes a density matrix as an input n[it is an iTensor MPS], an updates it 'in place' by the evolved one. Different Trotter orders are available. At order o=2 the error made at each time state is O(tau^3). At order o=3 the error made at each time state is O(tau^4.). At order o=4 the error made at each time state is O(tau^5).

* `io_util.h`: small and  simple methods for inputs and outputs (not specific to this type of simulations). It contains the class `Parameters`, which is used to handle a set of input parameters defined from a acommand line.

* `mps_mpo_utils.h` and `mps_mpo_utils.h`: basic but useful general methods for MPS and MPO (not specific to density matrices nor dissipative systems).
