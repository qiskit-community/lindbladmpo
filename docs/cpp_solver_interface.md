# The C++ solver interface  
This file describes the low-level C++ interface. Typical usage of the package does not require following most of the details described below.

One notable use-case is when running the solver executable on a remote server (or cluster node) with an input file generated using the Python interface on a different machine.

## C++ solver input files

The solver executable expects a single parameter in its command line,
which is the string literal `input_file`, followed by a space and the file name of an input file.
In the input file each parameter is specified on a separate line,
which contains the parameter name, a space, the equality sign '=',
a space, and a single value or a comma-separated list of values (without spaces).

The solver that is based on the ITesnor library, uses 1-based indexes for the qubits,
and this is reflected in the input and output files.
The Python code converts the 0-based indexes that are used by convention in Python indexing
when creating the input file and when loading the output files of the solver, so that the Python interface exposes 0-based indexes.

The parameters that the executable accepts have identical names and meaning
as the parameters detailed for the [Python interface](API_DOCS.md), with a few exceptions detailed in the following.  
* `b_unique_id`. This parameter is not recognized by the solver.
* `unique_id`, `metadata`. These string parameters are currently ignored by the solver, except for being allowed in the input file and saved in the log file.
* `J`, `J_z`. These two executable parameters are passed as a single scalar value
  (if the Python parameter is scalar), or as a list (no matrix format is supported).
  If either `J` or `J_z` is a matrix in the Python, then the list consists of
  all entries which are nonzero in either one of the two matrices.
  In this case, there are two additional parameters used in order to describe
  the qubit pairs to which the values correspond;
* `first_bond_indices`, `second_bond_indices`. Lists of the indices of the
  first and second qubit (respectively) of each entry in the lists
  `J`, `J_z` (that must be of identical length).

## C++ solver output files

The solver executable generates files: a saved state (at the simulation end), observables at requested times, global solver data at the same times as the observables, and a log file.
The format of the saved state file is not described here. 
The log file is textual. The observabl files have a consistent table format with each line holding one data value with tab-separated fields, as detailed in the following. We note that similarly to the input file, the output files qubit indexes are 1-based.

The three output files are:
* The one-qubit observables file, with the following structure as defined by the columns:
  `time`, `operator`, `index`, `value`. The time column stores the time at which the observable
  is taken, `operator` is the Pauli operator {X,Y,Z}, `index` is the qubit index and
  `value` is the Pauli expectation value.

* The two-qubit observables file has the columns:
  `time`, `operator`, `index_1`, `index_2`, `value`. 
  The `time` column stores the time at which the observable is taken,
  `operator` is the Pauli operator XX, XY, etc., `index_1` is the first qubit
  index, `index_2` is the second qubit index and `value` is the 
  two-qubit Pauli expectation value.

* The global quantities file has the columns: `time`, `quantity`, 
  `value`. The `time` column stores the time at which the quantity is taken,
  `quantity` is a string denoting the name of the calculated global quantity 
  and `value` is the value.

