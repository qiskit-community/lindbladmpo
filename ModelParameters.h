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
        operator[]("J") = "1";       //Hopping H=-J*(S+S- + S-S+) = -2*J*(SxSx+SySy). THis parameter can either be a single value, or a list of values
        operator[]("J_z") = "1";       //Sz-Sz Interaction strength. This parameter can either be a single value, or a list of values
        operator[]("h_x") = "0";     // magnetic field in the x direction. Note: h_x(sigma^+ + sigma^-)/2 = h_x*sigma^x/2 = h_x*S^x
        operator[]("h_y") = "0";     // magnetic field in the y direction. Note:                            h_y*sigma^y/2 = h_y*S^y
        operator[]("h_z") = "0";     // magnetic field in the z direction. Note:                            h_z*sigma^z/2 = h_z*S^z

        //Losses / dissipation
        operator[]("gamma") = "1.0"; //Strength of the "loss term"
        //Initial state
        operator[]("x_init") = "0"; //Initial state = spins pointing in the x direction
        operator[]("y_init") = "0"; //Initial state = spins pointing in the y direction

        //Lattice
        operator[]("b_periodic_x") = "false"; // if true -> periodic boundary conditions in the x direction (Warining: potential huge cost in terms of bond dimension)
        operator[]("b_periodic_y") = "false"; // if true -> periodic boundary conditions in the y direction
        operator[]("Lx") = "4";               //Small system by default
        operator[]("Ly") = "1";
        //Instead of the square lattice defined by tha above options, one can provide an explicit list of the couplings between q-bits
        //
        operator[]("N") = "0"; //Number of q-bits in the case of a used-defined lattice. If 0 => square lattice defined by the parameters above
        operator[]("Asites") = ""; //List of site indices
        operator[]("Bsites") = ""; //List of site indices


        //1-Qbit observables
        operator[]("1q_components") = "x,y,z"; // Vector of components
        operator[]("1q_sites") = "";           // Vector of long integers. If left empty => equivalent to 1,2,3,...,N

        //2-Qbit observables
        // Vector of components. Each element should have two letters. For instance XY means that <sigma^x(i)sigma^y(j)> will be computed for the pairs i,j specified in the argument "2q_sites".
        operator[]("2q_components") = "xx,yy,zz";
        operator[]("2q_sites") = ""; // Vector of long integers i1,j1,i2,j2,.... If left empty => equivalent to all pairs 1,2,1,3,...,1,N,    2,1,2,3,2,4,...,2,N,  ...  N,N-1
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
};

//____________________________________________________________________
#endif