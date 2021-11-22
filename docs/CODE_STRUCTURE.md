# A few words about the code structure

* The files `Pauli.cc` and `Pauli.h` contain the implementation of the vectorization for spin-1/2 systems, that is the representation of the many-body density matrix using an MPS. It is specific to 2-level systems (qubits), but not specific the special form of the Lindbadian of the specific model to be studied. It makes use of a few low-level routines of the iTensor library. Remark: The identity matrix plus the three Pauli matrices form a convenient basis of the the space of density matricxes of one qubit, hence the name of these two files.

* `lindbladmpo.cc`: contains the main() function. This is where the initial state is constructed (the way the initial state is defined depends on the input parameters), and where the loop over the time steps is defined. It is also where the observable of interests are computed, and the input and outputs (simple prints to the standard output, or to files) are handled.

Possible initial states: (1) A pure (product) state where all the spins point in a well-defined direction, or (2) An inifinite-temperature (mixed) state, where the density matrix is proportionnal to the identity matrix [to be supported soon], and (3) A mixed state which is read from the disk, after being saved at the end of some previous simulation.

* `lindbladian.h`: This is the place where the Lindbladian super-operator of the specific model to be simulated is defined. It takes the form of an `autoMPO` object of the iTensor library. Thanks to the use of the iTensor library, the terms in the Lindbladian can be defined in a simple way (using Pauli operators and products of such operators), almost as if one were writing them with pen and paper.

* `ModelParameters.h`: defines the parameters which are specific to the model to be simulated. These parameters are couples of the type "name" / "value" (values are floating point numbers, integers or strings). Their values can be defined in the command-line when calling the executable. Example `./lindbladmpo param1 val1 param2 val2`. Each parameter not specified in the command line will take its default value. The default values are defined in this file as well. The same arguments can also be specified using an input file. Example:  `./lindbladmpo input_file in.txt` where the file `in.txt` contains separate lines of the form `param1 = val1` etc.

* `SimulationParameters.h`: defines the parameters that are specific to the MPS/MPO simulation method, and in particular about the wanted level of  accuracy (maximum truncation, maximum bond dimension, etc.)

* `SimpleSquareLattice.h`: contains a basic class to encode the spatial geometry of the system (a lattice). In the present version it can handle a simple one-dimensional chain, or a square lattice with cylindrical boundary conditions. Similar classes could be created to handle other geometries and other lattices.

* `TimeEvolution.h` and `TimeEvolution.cc`: Contain the `TimeEvolver` class. Such class stores the parameters associated to a single-time-step evolution of the density matrix. An important parameter is for instance the value/length `tau` of one  time step. A `TimeEvolver` also contains other parameters associated to the approximations (truncations, etc.) to be made when applying such the time evolution operator (which is an MPO) to a given density matrix. The `TimeEvolver` class is independent of the details of the specific model to be studied. The actual time-evolution is coded in `TimeEvolution.cc`: the method `evolve` takes a density matrix as an input [it is an iTensor MPS], an updates it 'in place' by the evolved one. Different Trotter orders are available. At order o=2 the error made at each time state is O(tau^3). At order o=3 the error made at each time state is O(tau^4.). At order o=4 the error made at each time state is O(tau^5).

* `io_util.h`: small and  simple methods for inputs and outputs (not specific to this type of simulations). It contains the class `Parameters`, which is used to handle a set of input parameters (each one being a pair "name" / "value") defined from a acommand line.

* `mps_mpo_utils.h` and `mps_mpo_utils.h`: basic but useful general methods for MPS and MPO (not specific to density matrices nor dissipative systems).