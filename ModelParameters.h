#ifndef _MODELPARAMETERS_
#define _MODELPARAMETERS_

#include "SimulationParameters.h"
using namespace itensor;
using namespace std;

//____________________________________________________________________
class ModelParameters : public SimulationParameters
{
public:
    ModelParameters() : SimulationParameters() //We first call the SimulationParameters constructor
    {
        //Specify below all the allowed parameter names,
        //and their default values

        //Parameters of the Hamiltonian part of the model [Same notations & conventions as in https://doi.org/10.1103/PhysRevA.93.023821
        operator[]("U") = "0";       //Sz-Sz Interaction strength U*Sz*Sz
        operator[]("J") = "1";       //Hopping H=-J*(S+S- + S-S+) = -2*J*(SxSx+SySy)
        operator[]("Omega") = "0.5"; //magnetic field in the x direction. Note: Omega(sigma^+ + sigma^-) = Omega*sigma^x = 2*Omega*S^x
        operator[]("Delta") = "0";   // magnetic field in the z direction  H=Delta*S^z

        //Losses / dissipation
        operator[]("gamma") = "1.0"; //Strength of the "loss term"
        //Initial state
        operator[]("x_init") = "0"; //Initial state = spins pointing in the x direction
        operator[]("y_init") = "0"; //Initial state = spins pointing in the y direction

        //Lattice
        operator[]("b_periodic_x") = "false"; // if true -> periodic boundary conditions in the x direction (Warining: potential huge cost in terms of bond dimension)
        operator[]("b_periodic_y") = "false"; // if true -> periodic boundary conditions in the y direction
        operator[]("Lx") = "4";
        operator[]("Ly") = "1";

        operator[]("outputfilename") = "auto"; //if "auto", a name is generated using the most inportant parameters of the model (see FileName() below)
        
    }
    void check()
    {
        int Lx = longval("Lx");
        int Ly = longval("Ly");
        const int N = Lx * Ly;

        if (Lx < 1 || Ly < 1)
            cerr << "Error in the lattice parameters Lx=" << Lx << " Ly=" << Ly << endl, exit(1);
        if (N < 1)
            cerr << "Error, this code assumes that the system has at least 1 site.\n", exit(1);
    }
    string FileName() //This string can be used to construct the prefix of the output file names, so that it contain the most important parameters of the simulation
    {
        if (stringval("outputfilename") == "auto")
        {
            int Lx = longval("Lx");
            int Ly = longval("Ly");
            //double U = val("U");
            double J = val("J");
            double Omega = val("Omega");
            double Delta = val("Delta");
            double gamma = val("gamma");
            double tau = val("tau");
            string name = "Lx=" + to_string(Lx) + "_Ly=" + to_string(Ly) + "_Delta=" + to_string(Delta) + "_Omega=" + to_string(Omega) + "_J=" + to_string(J) + "_Gamma=" + to_string(gamma) //
                          + "_tau=" + to_string(tau) + "_chi=" + to_string(longval("MaxDimRho")) + "_Time=" + to_string(val("T"));

            if (longval("up_init") != 0)
                name += "_init=up";
            if (longval("down_init") != 0)
                name += "_init=down";
            if (longval("rho_inf_init") != 0)
                name += "_init=rho_inf";
            if (longval("x_init") != 0)
                name += "_init=x";
            if (longval("y_init") != 0)
                name += "_init=y";
            if (stringval("load_state_file") != "")
                name += "_init=fromfile";

            name += ".dat";
            return name;
        }
        else
        {
            return stringval("outputfilename");
        }
    }
};

//____________________________________________________________________
#endif