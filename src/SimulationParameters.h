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
        //The parameters below concern the simulation algorithm and approximations, these parameters are not model-specific
        //THe model-specific parametrs should be defined in ModelParameters.h

        operator[]("t_init") = "0";            //Initial time
        operator[]("t_final") = "1";            //Final time
        operator[]("tau") = "0.1";        //time step for the time evolution (should be small for better accuracy)
        operator[]("output_step") = "1";          //How often (in units of tau) shall we compute (and write to disk) the observables (just <Sz_i> etc. ). Set to zero if you do not want any observable to be computed

        operator[]("Trotter_order") = "4"; //Possible choices are 2,3,4. 3 or 4 are recommended.
        operator[]("max_dim_rho") = "400";  //maximum bond dimension for density matrices
        operator[]("min_dim_rho") = "1";  //minimum bond dimension for density matrices
        operator[]("cut_off_rho") = "1e-16"; //maximum truncation error for density matrices
        operator[]("b_force_rho_trace") = "1";  // We do rho/=trace{rho} at every time step, to keep trace[rho]=1, irrespectively of finite-step errors
        operator[]("b_force_rho_Hermitian") = "1"; // Replace rho by 0.5*(rho+rho^dagger) before measuring observables. This is recommanded, since it reduces some errors
        operator[]("b_initial_rho_orthogonalization") = "0";  //if <>0 = > After reading Rho from the disk, perform some orthogonalization/truncation

        operator[]("init_Pauli_state") = "+z";                                                                  
        //If not equal to "" the initial state (density matrix rho) will be read from the disk.
        //Give a filename if you want to read a previous density matrix from disk, as the initial state of the time evolution. 
        //Note that initial state is in fact stored in three files, with names *.state.ops, *.state.sites and *.state.rho
        operator[]("load_files_prefix") = "";
        //Whether the final state (density matrix) will be saved to disk.
        operator[]("b_save_final_state") = "0";   
                                              
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
