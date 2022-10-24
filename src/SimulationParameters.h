#ifndef _SIMULATIONPARAMETERS_
#define _SIMULATIONPARAMETERS_

using namespace itensor;
using namespace std;

//____________________________________________________________________
class SimulationParameters : public Parameters
{
public:
    SimulationParameters() //constructor
    {
        // The parameters below concern the simulation algorithm and approximations, and are not
        // model-specific (model-specific parameters should be defined in ModelParameters.h).

        operator[]("t_init") = "0";         // Initial time.
        operator[]("t_final") = "1";        // Final time.
        operator[]("tau") = "0.1";          // Step for the time evolution (should be smaller than
            // typical oscillation periods in the dynamics, but not too small for good performance)
        operator[]("output_step") = "1";    // Determines every how many tau time steps to compute
            // (and save) the observables. If set to 0, no observables are computed.

        operator[]("trotter_order") = "4";  // Possible choices are 2, 3, 4. 3 or 4 are recommended.
        operator[]("max_dim_rho") = "400";  // Maximum bond dimension for density matrices
        operator[]("cut_off_rho") = "1e-16"; // Maximum truncation error for density matrices. The actual
        	// truncation is done using the most severe condition between cut_off_rho and max_dim_rho.

        operator[]("b_force_rho_trace") = "1";		// Whether to force the density matrix trace to 1,
            // by substituting rho /= trace{rho} at every time step, compensating for finite-step errors
        operator[]("force_rho_hermitian_step") = "4";	// Determines every how many tau time steps
        	// to substitute rho = 0.5 * (rho + rho^dagger). This may reduce certain errors, but is
        	// computationally expensive.
        operator[]("b_initial_rho_compression") = "0";	// If nonzero, after reading rho from
         	// a saved file, perform a re-gauging/compression using ITensor's method orthogonalize().
        operator[]("b_quiet") = "0";	// If nonzero, after initialization of the simulation
         	// avoid the console output at every time step (but write it the log file).

        operator[]("init_product_state") = "";  // Initialize a product state. Can be specified
        	// for every qubit separately (as a comma-separated list), or uniformly for all qubits.
        	// The default, if left empty, is the "+z" Pauli state for all qubits, unless
        	// init_graph_state is nonempty, in which case all qubits are initialized to "+x".
        	// For a Pauli state, the format is + or - in the first character, and x, y, z in the
        	// second character, denoting a single-qubit state pointing along the positive/negative
        	// direction of the specified axis of the Bloch sphere.
        	// The string "id" is supported for indicating the fully mixed state.
        	// The format "p 0.25" is supported for indicating a diagonal density matrix, with
        	// p the probability of |0> (|0> = |up>).
        	// The format "q 0.7012 1.53" indicates the pure state:
        	// 0.7012 |0> + exp(i 1.53)sqrt{1-0.7012^2} |1>
        operator[]("init_pauli_state") = "";  // Initialize a Pauli state. This initialization
        	// is deprecated in favor of init_product_state, but remains currently supported
        	// as its synonym.
        operator[]("init_graph_state") = "";  // Initialize a graph state. A list of qubit
        	// pairs is expected, indicating all pairs to which a CZ gate is applied, after
        	// rotating all qubits (even if some are not included in the list) to the +x axis.
        	// Cannot be used together with any other initialization parameter.
        operator[]("init_cz_gates") = "";  // Initialize two-qubit CZ gates. A list of qubit
        	// pairs is expected, indicating all pairs to which a CZ gate is applied, after all
        	// qubits have been initialized arbitrarily according to the init_product_state
        	// parameter.

        operator[]("load_files_prefix") = "";	// If not an empty string, the initial state
         	// (density matrix rho) is to be read from the file system. Three files are being used,
         	// with names appended with ".state.ops", ".state.sites" and ".state.rho".
         	// This parameter cannot be used together with any other initialization parameter.
        operator[]("b_save_final_state") = "0";	// Whether to save the final state (density matrix) to
        	// the file system. Three files are generated, with the files names starting with the
        	// `output_files_prefix` string, with the endings ".state.rho", ".state.sites", and ".state.ops".
                                              
        operator[]("unique_id") = ""; // An optional unique id identifying the simulation. Not currently used (except for being saved in the input and log files).
        operator[]("metadata") = ""; // An optional user information, ignored by the solver (except for being saved in the input and log files).
        operator[]("input_file") = ""; // If not empty => name of the file from which some parameters must be read (in addition to the command line ones).
        operator[]("output_files_prefix") = "lindblad"; // Path and prefix of the file names where various simulation output is written

    }
};

        //OLD initial state settings
/*        operator[]("x_init") = "0";    //Initial state = spins pointing in the x direction
        operator[]("y_init") = "0";    //Initial state = spins pointing in the y direction
        operator[]("xm_init") = "0";    //Initial state = spins pointing in the -x direction
        operator[]("ym_init") = "0";    //Initial state = spins pointing in the -y direction
        operator[]("up_init") = "1";   //Set to "1" if you want to start from a (pure) state where all the spins are "up"
        operator[]("down_init") = "0"; //Set to "1" if you want to start from a (pure) state where all the spins are "up"
        operator[]("rho_inf_init") = "0";     //Set to "1" if you want to start from rho=identity (infinite temperature density matrix)

        operator[]("max_dim_psi") = "500";     //maximum bond dimension for wave-functions
        operator[]("cut_off_psi") = "1e-9";    //maximum truncation error  for wave-functions
        operator[]("sweeps") = "0";     //maximum number of sweeps in the DMRG (to get the initial wave-function as the ground-state of H0)
        operator[]("energy") = "1e-12"; //convergence criterium on the energy. This is only used in the initial DMRG runs to coinstruct (if needed) the initial state as the ground state of H0.

        operator[]("dmrgLDL") = "0";     //Set to "1" if you want to do some DMRG sweeps to minimize Ldagger*L (L:Lindbladian), to get (or approach) the steady state
        operator[]("LDLconv") = "5e-4";  //convergence criterium on L^dagger * L -- this is a relative value, so the DMRG stops when |dE/E|< LDLconv. In an exact calculation of the steady state E=0
        operator[]("sweepsRho") = "999"; //maximum number of sweeps in the DMRG (to get the initial density matrix as the ground-state of L^dagger*L)
        //If not equal to "" the initial state is a pure state (wavefunction) that will be read from the disk.
        operator[]("load_purestate_file") = "";
        //If not equal to "" the initial state is a pure state (wavefunction) that will be read from the disk.
        operator[]("save_purestate_file") = "";
*/
//____________________________________________________________________
#endif
