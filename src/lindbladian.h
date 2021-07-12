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

#ifndef _LINDBLADIAN_
#define _LINDBLADIAN_

using namespace itensor;
using namespace std;

//____________________________________________________________________

void SetLindbladian(SpinHalfSystem &C, ModelParameters param, Lattice2d L)
{
    // -----------------------------------------------------------
    // We first construct the Hamiltonian (unitary evolution) terms
    const unsigned int N = C.N;
    const unsigned int num_bonds = L.I.size();
    vector<double> h_x = param.doublevec("h_x");
    vector<double> h_y = param.doublevec("h_y");
    vector<double> h_z = param.doublevec("h_z");
    unsigned int h_x_len = h_x.size();
    unsigned int h_y_len = h_y.size();
    unsigned int h_z_len = h_z.size();

    if (h_x_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter h_x has " << h_x_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". h_x should be uniform for such a lattice.\n", exit(1);
    if (h_y_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter h_y has " << h_y_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". h_y should be uniform for such a lattice.\n", exit(1);
    if (h_z_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter h_z has " << h_z_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". h_z should be uniform for such a lattice.\n", exit(1);

    vector<double> g_0 = param.doublevec("g_0");
    vector<double> g_1 = param.doublevec("g_1");
    vector<double> g_2 = param.doublevec("g_2");
    unsigned int g_0_len = g_0.size();
    unsigned int g_1_len = g_1.size();
    unsigned int g_2_len = g_2.size();

    if (g_0_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter g_0 has " << g_0_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". g_0 should be uniform for such a lattice.\n", exit(1);
    if (g_1_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter g_1 has " << g_1_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". g_1 should be uniform for such a lattice.\n", exit(1);
    if (g_2_len != 1 && L.predefined && !L.predefined_chain)
        cerr << "Error: the paramter g_2 has " << g_2_len << " value(s) but L.predefined_chain=" << L.predefined_chain << ". g_2 should be uniform for such a lattice.\n", exit(1);

    if (h_x_len != 1 && h_x_len != N)
        cerr << "Error: the paramter h_x has " << h_x_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (h_y_len != 1 && h_y_len != N)
        cerr << "Error: the paramter h_y has " << h_y_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (h_z_len != 1 && h_z_len != N)
        cerr << "Error: the paramter h_z has " << h_z_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (h_x_len == 1)
        h_x = vector<double>(N, h_x[0]);
    if (h_y_len == 1)
        h_y = vector<double>(N, h_y[0]);
    if (h_z_len == 1)
        h_z = vector<double>(N, h_z[0]);

    vector<double> J = param.doublevec("J");
    vector<double> J_z = param.doublevec("J_z");

    if ((J_z.size() != 1 || J.size() != 1) && L.predefined)
    {
        cerr << "Error: J_z.size()=" << J_z.size() << " and " << J.size() << " but L.predefined=" << L.predefined << ". Couplings J and J_z must be uniform in such a predefined lattice.\n", exit(1);
    }

    if (J.size() != 1 && J.size() != num_bonds)
        cerr << "Error: the paramter J has " << J.size() << " values but 1 or " << num_bonds << " value(s) were expected.\n", exit(1);
    if (J.size() == 1)
    {
        J = vector<double>(num_bonds, J[0]);
    }

    if (J_z.size() != 1 && J_z.size() != num_bonds)
        cerr << "Error: the paramter J_z has " << J_z.size() << " values but 1 or " << num_bonds << " value(s) were expected.\n", exit(1);
    if (J_z.size() == 1)
    {
        J_z = vector<double>(num_bonds, J_z[0]);
    }

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
            auto_L += -.5 * J_z[n], "Sz", i, "Sz", j;

            auto_L += J[n], "_S+", i, "_S-", j;
            auto_L += J[n], "_S-", i, "_S+", j;
            auto_L += .5 * J_z[n], "_Sz", i, "_Sz", j;
        }

        //Magnetic field terms:
        for (int j = 1; j <= int(N); ++j)
        {
            if (h_x[j - 1] != 0.)
            {
                auto_L += .5 * h_x[j - 1], "Sx", j;
                auto_L += -.5 * h_x[j - 1], "_Sx", j;
            }
            if (h_y[j - 1] != 0.)
            {
                auto_L += .5 * h_y[j - 1], "Sy", j;
                auto_L += -.5 * h_y[j - 1], "_Sy", j;
            }
            if (h_z[j - 1] != 0.)
            {
                auto_L += .5 * h_z[j - 1], "Sz", j;
                auto_L += -.5 * h_z[j - 1], "_Sz", j;
            }
        }
    }
    // -----------------------------------------------------------
    // Dissipative terms

    if (g_0_len != 1 && g_0_len != N)
        cerr << "Error: the paramter g_0 has " << g_0_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (g_1_len != 1 && g_1_len != N)
        cerr << "Error: the paramter g_1 has " << g_1_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (g_2_len != 1 && g_2_len != N)
        cerr << "Error: the paramter g_2 has " << g_2_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
    if (g_0_len == 1)
        g_0 = vector<double>(N, g_0[0]);
    if (g_1_len == 1)
        g_1 = vector<double>(N, g_1[0]);
    if (g_2_len == 1)
        g_2 = vector<double>(N, g_2[0]);
    // OLD: const double gamma = param.val("gamma");
    // cout << "Strength of the Lindblad terms, gamma=" << gamma << endl;

    for (int i = 1; i <= int(N); i++)
        // AddSingleSpinBath:
        // The first argument is the rate of dissipative processes where a spin goes from down to up
        // The second argument is the rate of dissipative processes where a spin goes from up to down
        // The third argument is the rate of energy-conserving (pure) dephasing processes
        C.AddSingleSpinBath(g_0[i - 1], g_1[i - 1], g_2[i - 1], i);
}
//____________________________________________________________________

#endif
