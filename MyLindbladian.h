#ifndef _MYLINDBLADIAN_
#define _MYLINDBLADIAN_

using namespace itensor;
using namespace std;

//____________________________________________________________________

void SetLindbladian(SpinHalfSystem &C, ModelParameters param, Lattice2d L)
{
    // -----------------------------------------------------------
    // We first construct the Hamiltonian (unitary evolution) terms
    const unsigned int N = C.N;
    const unsigned int num_bonds=L.I.size();
    vector<double> h_x = param.doublevec("h_x");
    unsigned int h_x_len = h_x.size();
    if (h_x_len != 1 && h_x_len != N)
        cerr << "Error: the paramter h_x has " << h_x_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (h_x_len == N)
    {
        h_x = vector<double>(N, h_x[0]);
    }
    vector<double> J = param.doublevec("J");
    if (J.size() != 1 && J.size() != num_bonds)
        cerr << "Error: the paramter J has " << J.size() << " values  but  1 or " << num_bonds << " value(s) were expected.\n", exit(1);
    if (J.size() == 1)
    {
        J = vector<double>(num_bonds, J[0]);
    }
    vector<double> J_z = param.doublevec("J_z");
    if (J_z.size() != 1 && J_z.size() != num_bonds)
        cerr << "Error: the paramter J has " << J_z.size() << " values  but  1 or " << num_bonds << " value(s) were expected.\n", exit(1);
    if (J_z.size() == 1)
    {
        J_z = vector<double>(num_bonds, J_z[0]);
    }

    const double Delta = param.val("Delta");

    cout << "SetHamiltonian with N=" << N << endl;
    for (int dag = 0; dag <= 1; dag++)
    {
        AutoMPO &auto_L = (dag == 0) ? (C.Lindbladian) : (C.LindbladianDag); //Linbladian operator, or its hermitian conjugate (samething for the unitary terms)

        //Note about the Lindblad equation: since we evolve rho (and not a wave function), each term in H
        //acts once with a "+" on the right of rho, and once with a "-" to the left of rho (prefix "_" in the operator name)

        for (unsigned int n = 0; n < num_bonds; n++)
        {
            int i = L.I[n], j = L.J[n];
            auto_L += -J[n], "S+", i, "S-", j;
            auto_L += -J[n], "S-", i, "S+", j;
            auto_L += J_z[n], "Sz", i, "Sz", j;

            auto_L += J[n], "_S+", i, "_S-", j;
            auto_L += J[n], "_S-", i, "_S+", j;
            auto_L += -J_z[n], "_Sz", i, "_Sz", j;
        }

        //Magnetic field terms:
        for ( int j = 1; j <= N; ++j)
        {
            auto_L += h_x[j - 1], "Sx", j;
            auto_L += Delta, "Sz", j;

            auto_L += - h_x[j - 1], "_Sx", j;
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
