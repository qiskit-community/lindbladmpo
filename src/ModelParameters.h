// Licensed under the Apache License, Version 2.0 (the "License"); you may
// not use this file except in compliance with the License. You may obtain
// a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.

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
        operator[]("N") = "4"; // Number of qubits in the case of a used-defined lattice.
        // If N = 0, the total number of qubits is taken to be l_x * l_y (which must be both nonzero).
        // The Python interface does not provide a default value (requires a user assignment), but here
        // we assign a default value of 4 to allow running the solver from the command line for testing.

        //Parameters of the single-qubit Hamiltonian part of the model
        operator[]("h_x") = "0"; // magnetic field in the x direction. Note: h_x(sigma^+ + sigma^-)/2 = h_x*sigma^x/2 = h_x*S^x
        operator[]("h_y") = "0"; // magnetic field in the y direction. Note:                            h_y*sigma^y/2 = h_y*S^y
        operator[]("h_z") = "0"; // magnetic field in the z direction. Note:                            h_z*sigma^z/2 = h_z*S^z

        //Losses / dissipation
        operator[]("g_0") = "0"; // Strength of the excitation term
        operator[]("g_1") = "0"; // Strength of the loss term
        operator[]("g_2") = "0"; // Strength of the dephasing term (a Pauli z)
        operator[]("g_3") = "0"; // Strength of the bit flip term (a Pauli x)
        operator[]("g_4") = "0"; // Strength of the bit-phase flip term (a Pauli y)

        //Parameters of the interaction Hamiltonian
        operator[]("J") = "0";   // Hopping H=-J*(S+S- + S-S+) = -2*J*(SxSx+SySy). THis parameter can either be a single value, or a list of values
        operator[]("J_z") = "0"; // Sz-Sz Interaction strength. This parameter can either be a single value, or a list of values

        //Lattice specification
        operator[]("b_periodic_x") = "false"; // if true -> periodic boundary conditions in the x direction (Warining: potential huge cost in terms of bond dimension)
        operator[]("b_periodic_y") = "false"; // if true -> periodic boundary conditions in the y direction
        operator[]("l_x") = "0"; // The default setting of l_x = 0, l_y = 1 will construct a 1D chain of length N.
        operator[]("l_y") = "1";
        // Instead of the lattice defined by tha above options, the user can provide an explicit list of the couplings between qubits
        operator[]("first_bond_indices") = ""; // List of indices, first site of each bond
        operator[]("second_bond_indices") = ""; // List of indices, second site of each bond

        //1-qubit observables
        operator[]("1q_components") = "z"; // Vector of components
        operator[]("1q_indices") = ""; // Vector of integers. If left empty => equivalent to 1,2,3,...,N

        //2-qubit observables
        // Vector of components. Each element should have two letters. For instance XY means that <sigma^x(i)sigma^y(j)> will be computed for the pairs i,j specified in the argument "2q_indices".
        operator[]("2q_components") = "zz";
        operator[]("2q_indices") = ""; // Vector of integers i1,j1,i2,j2,.... If left empty => equivalent to all pairs 1,2,1,3,...,1,N,    2,1,2,3,2,4,...,2,N,  ...  N,N-1

        //3-qubit observables
        // Vector of components. Each element should have three letters.
        operator[]("3q_components") = "";
        operator[]("3q_indices") = ""; // Vector of integers i1,j1,k1,i2,j2,k2,.... If left empty nothing is calculated!

        operator[]("custom_observables") = "";
        operator[]("collapse") = "";
    }
};

//____________________________________________________________________
#endif