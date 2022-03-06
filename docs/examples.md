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
  


