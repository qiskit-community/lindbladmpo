#ifndef _MODELPARAMETERS_
#define _MODELPARAMETERS_

#include "utils.h"
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
        operator[]("N") = 6;        //Number of spins
        operator[]("periodic") = 0; // if not zero -> periodic boundary conditions in the x direction (Warining: potential huge cost in terms of bond dimension)
        operator[]("U") = 0;        //Sz-Sz Interaction strength U*Sz*Sz
        operator[]("J") = 1;        //Hopping H=-J*(S+S- + S-S+) = -2*J*(SxSx+SySy)
        operator[]("Omega") = 0.5;  //magnetic field in the x direction. Note: Omega(sigma^+ + sigma^-) = Omega*sigma^x = 2*Omega*S^x
        operator[]("Delta") = 0;    // magnetic field in the z direction  H=Delta*S^z

        //Losses / dissipation
        operator[]("gamma") = 1.0; //Strength of the "loss term"
        //Initial state
        operator[]("x_init") = 0; //Initial state = spins pointing in the x direction
        operator[]("y_init") = 0; //Initial state = spins pointing in the y direction

        //2d cylinders:
        operator[]("cylinder") = 0; // If 0 => chain (defined by N only); if <>0 => cylinder of dimension Lx * Ly =N
        operator[]("Lx") = 0;
        operator[]("Ly") = 0;
    }
    void check()
    {
        const int N = longval("N");
        int Lx = longval("Lx");
        int Ly = longval("Ly");
        if (longval("cylinder"))
        {
            if (Lx < 1 || Ly < 1)
                cout << "Error: cylinder=" << longval("cylinder") << " Lx=" << Lx << " Ly=" << Ly << endl, exit(0);
            operator[]("N") = Lx * Ly;
        }
        else
        {
            if (Lx != 0 || Ly != 0)
                cout << "Error: cylinder=" << longval("cylinder") << " Lx=" << Lx << " Ly=" << Ly << endl, exit(0);
            operator[]("Ly") = 1;
            operator[]("Lx") = N;
        }
        if (N < 1)
            cout << "Error, this code assumes that the system has at least 1 site.\n", exit(0);
    }
    string FileName() //This string can be used to construct the output data file name, so that it contain the most important parameters of the simulation
    {
        int N = longval("N");
        int Lx = longval("Lx");
        int Ly = longval("Ly");
        double U = val("U");
        double J = val("J");
        double Omega = val("Omega");
        double Delta = val("Delta");
        double gamma = val("gamma");
        double tau = val("tau");
        return "_N=" + to_string(N) + ((Lx > 0) ? ("_Lx=" + to_string(Lx)) : ("")) + ((Ly > 0) ? ("_Ly=" + to_string(Ly)) : (""))  //
               + "_Delta=" + to_string(Delta) + "_Omega=" + to_string(Omega) + "_J=" + to_string(J) + "_Gamma=" + to_string(gamma) //
               + "_tau=" + to_string(tau) + "_chi=" + to_string(longval("MaxDimRho")) + "_Time=" + to_string(val("T"));
    }
};

//____________________________________________________________________
#endif