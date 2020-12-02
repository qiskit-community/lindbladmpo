#ifndef _MYLINDBLADIAN_
#define _MYLINDBLADIAN_

#include "utils.h"
using namespace itensor;
using namespace std;

//____________________________________________________________________

void SetLindbladian(SpinHalfSystem &C, ModelParameters param, Lattice2d L)
{
    // -----------------------------------------------------------
    // We first construct the Hamiltonian (unitary evolution) terms
    const int N = C.N;
    const double U = param.val("U");
    const double J = param.val("J");
    const double Omega = param.val("Omega");
    const double Delta = param.val("Delta");

    cout << "SetHamiltonian: N=" << N << "\t Delta=" << Delta << "\tU=" << U << " J=" << J << " Omega=" << Omega << endl;
    for (int dag = 0; dag <= 1; dag++)
    {
        AutoMPO &auto_L = (dag == 0) ? (C.Lindbladian) : (C.LindbladianDag); //Linbladian operator, or its hermitian conjugate (samething for the unitary terms)

        //Note about the Lindblad equation: since we evolve rho (and not a wave function), each term in H
        //acts once with a "+" on the right of rho, and once with a "-" to the left of rho (prefix "_" in the operator name)

        for (unsigned int n = 0; n < L.I.size(); n++)
        {
            int i = L.I[n], j = L.J[n];
            auto_L += -J, "S+", i, "S-", j;
            auto_L += -J, "S-", i, "S+", j;
            auto_L += U, "Sz", i, "Sz", j;

            auto_L += J, "_S+", i, "_S-", j;
            auto_L += J, "_S-", i, "_S+", j;
            auto_L += -U, "_Sz", i, "_Sz", j;
        }

        //Uniform magnetic field everywhere:
        for (int j = 1; j <= N; ++j)
        {
            auto_L += 2 * Omega, "Sx", j;
            auto_L += Delta, "Sz", j;

            auto_L += -2 * Omega, "_Sx", j;
            auto_L += -Delta, "_Sz", j;
        }
    }
    // -----------------------------------------------------------
    // Dissipative terms

    const double gamma = param.val("gamma");
    cout << "Strength of the Lindblad terms, gamma=" << gamma << endl;

    for (int i = 1; i <= N; i++)
        // AddSingleSpinBath:
        // The first argument (here =0) is related to dissipative processes where a spin goes from down to up
        // The 2n  argument (here = gamma * 0.5) is related to dissipative processes where a spin goes from up  to down
        C.AddSingleSpinBath(0, gamma * 0.5, i);
}
//____________________________________________________________________

#endif
