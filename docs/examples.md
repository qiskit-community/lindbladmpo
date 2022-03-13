# Simulation examples

The `/lindbladmpo/examples` directory contains three subdirectories with Python example files that can be useful for users interested in performing research using the MPO solver. 

* `/examples/disordered_chain` - Contains a single file, [disordered_chain.py](/lindbladmpo/examples/disordered_chain/disordered_chain.py),
  with a basic example that is relatively similar to the tutorial in 
  [qubit_chain_tutorial.ipynb](qubit_chain_tutorial.ipynb), which can be
  used as another starting point for using the solver. The problem being
  set in this example is of a chain of qubits with the middle qubit being
  driven, and according to a disorder parameter, resonant spreading of induced 
  excitations or their localization can be observed.
    
* `/examples/simulation_building` - The main class in this example is
  [`LindbladMatrixSolver`](/lindbladmpo/examples/simulation_building/LindbladMatrixSolver.py)
  that implements an interface almost identical to
  that of the MPO solver class (`LindbladMPOSolver`, from which it is derived), and employs 
  [qiskit-dynamics](https://github.com/Qiskit/qiskit-dynamics)
  to solve the same dynamical problem using a scipy ode integrator. Classes
  facilitate writing Hamiltonian and Lindbladian operators in code are
  being used for the creation of `numpy` matrices for the solver.
  This example is relied upon in the third, most extensive example;
  
* `/examples/qubit_driving` - Contains an example for a large scale
  research project performed using the solver package. The entry point file (to run) is 
  [run-driving.py](/lindbladmpo/examples/qubit_driving/run-driving.py),
  which runs a simulation according some parameters given in this file,
  saves a dataframe entry describing the simulation using its metadata
  (that can be later queried when there are multiple simulation), and saves 
  the output data and figure files when the simulation returns. The simulation
  can be easily run using the MPO solver and also using a brute-force (full
  density matrix) simulation using [qiskit-dynamics](https://github.com/Qiskit/qiskit-dynamics)
  (for low enough qubit numbers),
  to compare and validate the numerical results.
  
In the disordered chain and the qubit driving project examples above, there is an extensive usage of plotting functions imported from the file [`plot_routines.py`](/lindbladmpo/plot_routines.py), which are very useful for easily plotting the solver output (observables and entropies) in different formats - in curves as a function of time, connected correlation matrices, and space-time (qubit-time) diagrams.

## Generating an input file for simulations on a remote machine
  
An important use case is to run the solver executable on a remote server (or cluster node) with an input file generated using the Python interface on a different machine.  In this case, follow the steps below after generaring an inout file using the Python interface:
* Update the paths inserted in the input file to be those relevant for target machine.
* Run the solver from the command line on the target machine. The solver executable expects a single parameter in its command line,
which is the string literal `input_file`, followed by a space and the file name of the input file. See the installation instructions for guaranteeing the proper multithreaded BLAS/LAPACK libraries are available.
* When the solver terminates, use the Python interface to load the output data files and plot the results.


