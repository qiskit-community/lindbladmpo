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

#include "itensor/all.h"
#include "Pauli.h"
#include "io_utils.h"
#include "mps_mpo_utils.h"
#include "SimpleSquareLattice.h"
#include "TimeEvolution.h"
#include "ModelParameters.h"
#include "lindbladian.h"
#include "gates.h"
#include <chrono>
#include <sstream>
#include <iostream>

using namespace itensor;
using namespace std;
using namespace std::chrono;

stream2d cout2 = stream2d(&cerr, NULL);

const string SOLVER_VERSION = "0.2.4";

const double IMAGINARY_THRESHOLD = 1e-4;
// Threshold for the imaginary value of a quantity that should be real, to issue a warning

const double TRACE_RHO_THRESHOLD = 1e-4;
// Threshold for the deviation of the density matrix trace from 1, to issue a warning


void validate_2q_list(vector<long> &vect, int N, string const &list_name);
void validate_3q_list(vector<long> &vect, int N, string const &list_name);

int main(int argc, char *argv[])
{
	ModelParameters param;

	auto t_start_sim = steady_clock::now();
	// Start tracking simulation time

	// Read the input parameters given in the command line. Default values are substituted for all
	// missing parameters, as defined in the files SimulationParameters.h and ModelParameters.h
	param.ReadArguments(argc, argv);

	// If a filename was given as the parameter `input_file`, read the parameters given in this file.
	string inputfilename = param.stringval("input_file");
	if (inputfilename != "")
		param.ReadFromFile(inputfilename);
	// Now `param` contains all parameters.

	Lattice2d lattice;
	int N = param.longval("N");
	int Lx = param.longval("l_x");
	int Ly = param.longval("l_y");
	vector<long> A = param.longvec("first_bond_indices");
	vector<long> B = param.longvec("second_bond_indices");
    const unsigned int n_bonds = A.size();
    const double _2_N = pow(2, -(double)N);

	stringstream strstr = stringstream();
	if (n_bonds)
	{
		// A user-defined lattice given explicitly using the bond couplings. In this case,
		// the number of qubits is taken to be N (validated in ``Lattice2d``, and l_x, l_y
		// parameters must take their default values, otherwise it's inconsistent
		if (!(Lx == 0 && Ly == 1))
			cerr << "Error, when bond couplings are specified explicitly, the parameters "
				"(l_x, l_y) must be left at their default values (0, 1).\n", exit(1);
		lattice = Lattice2d(N, A, B, strstr);
	}
	else
	{
		// A square lattice
		if (Lx > 0 && Ly > 0)
		{
			// If N is nonzero, it has to be consistent with l_x, l_y
			if (N > 0 && N != Lx * Ly)
				cerr << "Error, invalid N = " << N << ", not equal to l_x * l_y.\n", exit(1);
			N = Lx * Ly; // If N is zero, it is assigned with l_x * l_y
		}
		else if (Lx == 0 && Ly == 1)
			Lx = N; // The Lattice2d() constructor below will verify N.
		else
			cerr << "Error, invalid (l_x, l_y) = " << Lx << "," << Ly << ".\n", exit(1);
		lattice = Lattice2d(Lx, Ly, strstr, param.boolval("b_periodic_x"), param.boolval("b_periodic_y"));
	}
	string output_prefix = param.stringval("output_files_prefix");
	ofstream log_file(output_prefix + ".log.txt");
	cout2 = stream2d(&cout, &log_file);

	cout2.precision(8);
	cout2 << "lindbladmpo solver log. Solver version: " << SOLVER_VERSION << "\n";
    cout2 << "---------------------------------------------\n";
	param.Print(cout2);
    cout2 << "------------------------------------------------------------------\n";
	cout2 << strstr.str();
    cout2 << "------------------------------------------------------------------\n";

	SpinHalfSystem C(N);

	string load_prefix = param.stringval("load_files_prefix");
	if (load_prefix != "")  // || param.stringval("load_purestate_file") != "")
	{
//      if (param.stringval("load_purestate_file") != "")
//          file_name = param.stringval("load_purestate_file");
		string file_name = load_prefix;
		string f1 = file_name + ".state.ops";
		string f3 = file_name + ".state.sites";
		cout2 << "Opening '" << f1 << "' and '" << f3 << "'...";
		cout2.flush();
		readFromFile(f3, C.sites);
		readFromFile(f1, C.siteops);
		cout2 << "done.\n";
		C.rho = MPS(C.siteops);
		// Rho is still undefined but its structure is initialized using C.siteops
		C.Lindbladian = AutoMPO(C.siteops);
	}

    C.ConstructIdentity(); // Construct the density matrix corresponding to infinite temperature

	// Note on the option below: if we  "Normalize=true" (default for fitapplyMPO)
	// we would normalize rho such that Tr[rho^2]=1, (norm of the MPS)
	// which is of course not appropriate for a density matrix.
	auto argsRho = Args();
	argsRho.add("Normalize", false);
	argsRho.add("MaxDim", param.longval("max_dim_rho"));
	argsRho.add("Cutoff", param.val("cut_off_rho"));

	vector<long> init_graph_state = param.longvec("init_graph_state");
	vector<long> init_cz_gates = param.longvec("init_cz_gates");
	bool b_graph_state = false, b_cz_pairs = false;
	string s_cz_param = "";
	if (init_graph_state.size() > 0)
	{
		b_cz_pairs = true;
		b_graph_state = true;
		s_cz_param = "init_graph_state";
		if (init_cz_gates.size() == 0)
			init_cz_gates = init_graph_state;
		else
			cout2 << "Error: the parameter init_cz_gates cannot be used " <<
				"if init_graph_state is nonempty.\n", exit(1);
	}
	else if (init_cz_gates.size() > 0)
	{
		b_cz_pairs = true;
		s_cz_param = "init_cz_gates";
	}
	if (init_cz_gates.size() % 2 == 1)
		cout2 << "Error: the list of indices given in the parameter " << s_cz_param <<
			" should have an even length.\n", exit(1);

	vector<string> init_pauli_state = param.stringvec("init_pauli_state");
	vector<string> init_product_state = param.stringvec("init_product_state");
	string s_init_param;
	if (init_pauli_state.size() != 0)
	{
		s_init_param = "init_pauli_state";
		cout2 << "Warning: the init_pauli_state parameter has been deprecated and " <<
		 	"will be removed in the future. Please use init_product_state instead.\n";
		if (init_product_state.size() == 0)
			init_product_state = init_pauli_state;
		else
			cout2 << "Error: the parameter init_pauli_state cannot be used " <<
				"if init_product_state is nonempty.\n", exit(1);
	}
	else
		s_init_param = "init_product_state";
	unsigned int init_len = init_product_state.size();

	if (load_prefix != "" && (b_cz_pairs || init_len > 0))
		cout2 << "Error: if load_files_prefix is nonempty, no other initialization " <<
			"parameter can be used.\n", exit(1);
	if (init_len == 0)
	{
		if (b_graph_state)
			init_product_state = vector<string>(N, "+x");
		else
			init_product_state = vector<string>(N, "+z");
	}
	else
	{
		if (b_graph_state)
			cout2 << "Error: If init_graph_state is nonempty, no other initialization " <<
				"parameter can be used.\n", exit(1);
		if (init_len == 1)
			init_product_state = vector<string>(N, init_product_state[0]);
		else if (int(init_len) != N)
			cout2 << "Error: the parameter " << s_init_param << " has " << init_len <<
				" value(s) but 0, 1 or " << N << " value(s) were expected.\n", exit(1);
	}
	if (b_cz_pairs)
		validate_2q_list(init_cz_gates, N, s_cz_param);

	const double t_0 = param.val("t_init");
	const double tau = param.val("tau");
	const double t_f = param.val("t_final");
	const double t_total = t_f - t_0;
	MPS psi;
	bool psi_defined = false;
	if (load_prefix != "")
	{ //Read the density matrix from disk
	    string file_name = load_prefix;
		file_name += ".state.rho";
		cout2 << "Read the initial rho from the file '" << file_name << "'...";
		cout2.flush();
		readFromFile(file_name, C.rho);
		cout2 << "done.\n";
		cout2 << "Bond dimension of rho:" << maxLinkDim(C.rho) << "\n";
	}
	else
	{
    	// Set the initial wavefunction matrix product state
		auto initState = InitState(C.sites);
		for (int i = 1; i <= N; ++i) // Start with all spins up
			initState.set(i, "Up");
		psi = MPS(initState);
		vector<double> a_mixed_state = vector<double>(N, -1.);
		vector<double> dr_mixed_state = vector<double>(N, 0.);
		vector<double> di_mixed_state = vector<double>(N, 0.);

		// We init values outside of the valid range, that will be ignored below unless overwritten
		double sqrt05 = pow(.5, .5);
		for (int site_number = 1; site_number <= N; site_number++)
        {
        	const unsigned int i_site = (unsigned int)site_number - 1;
			string s_init = init_product_state[i_site];
			double R0 = .0, I0 = .0, R1 = .0, I1 = .0;
			if (s_init == "+z")
				R0 = 1.;
			else if (s_init == "-z")
				R1 = 1.;
			else if (s_init == "+x")
				R0 = R1 = sqrt05;
			else if (s_init == "+y")
				R0 = I1 = sqrt05;
			else if (s_init == "-x")
			{
			    R0 = sqrt05;
				R1 = -sqrt05;
			}
			else if (s_init == "-y")
			{
			    R0 = sqrt05;
				I1 = -sqrt05;
			}
			else
			{
				if (s_init == "id")
				{
					R0 = 1.;  // A valid value must be set to get a regular pure state
					a_mixed_state[i_site] = .5;
				}
				else
				{
					string s_type = "";
					if (s_init.size() > 2)
					 	s_type = s_init.substr(0, 2);
					if (s_type == "p ")
					{
						double b = 1.;  // Value is not expected to be used.
						string s_val = s_init.substr(2);
						try
						{
							b = stod(s_val);
						}
						catch (...)
						{
							cout2 << "Error: " << s_init << " is an unknown 1-qubit initial state. "
								"After the prefix \"p \" should follow a valid float value.\n", exit(1);
						}
						if (b < 0. || b > 1.)
							cout2 << "Error: " << s_val << " is an unknown 1-qubit initial "
								"mixed state coefficient (should be in the range [0, 1]).\n", exit(1);
						R0 = 1.;  // A valid value must be set to get a regular pure state
						a_mixed_state[i_site] = b;
					}
					else if (s_type == "q ")
					{
						string s_val;
						float val = 0.;
						vector<double> values;
						bool b_invalid = false;
						istringstream s_stream(s_init.substr(2));
						while (getline(s_stream, s_val, ' '))
						{
							try
							{
								val = stod(s_val);
							}
							catch (...)
							{
								b_invalid = true;
								break;
							}
							values.push_back(val);
						}
						if (values.size() != 2)
							b_invalid = true;
						if (b_invalid)
							cout2 << "Error: " << s_init << " is an unknown 1-qubit initial state. "
								"After the prefix \"q \" should follow two space-separated floats.\n",
								exit(1);
						R0 = cos(values[0] / 2.);
						double r1 = pow(1. - pow(R0, 2), .5);
						R1 = r1 * cos(values[1]);
						I1 = r1 * sin(values[1]);
					}
					else if (s_type == "r ")
					{
						string s_val;
						float val = 0.;
						vector<double> values;
						bool b_invalid = false;
						istringstream s_stream(s_init.substr(2));
						while (getline(s_stream, s_val, ' '))
						{
							try
							{
								val = stod(s_val);
							}
							catch (...)
							{
								b_invalid = true;
								break;
							}
							values.push_back(val);
						}
						if (values.size() != 3)
							b_invalid = true;
						if (b_invalid)
							cout2 << "Error: " << s_init << " is an unknown 1-qubit initial state. "
								"After the prefix \"r \" should follow three space-separated floats.\n",
								exit(1);
						R0 = 1.;
						a_mixed_state[i_site]=values[0];
						dr_mixed_state[i_site]=values[1];
						di_mixed_state[i_site]=values[2];
					}
					else
					{
						cout2 << "Error: " << s_init << " is an unknown 1-qubit initial state. "
							"It should be a string in {+x, -x, +y, -y, +z, -z, id}, \"p \" "
							"followed by a float, or \"q \" followed by two space-separated floats).\n",
							exit(1);
					}
				}
			}

            auto spin_ind = siteIndex(psi, site_number);
            if (site_number == 1)
			{
				Index ri = rightLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, ri = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, ri = 1, Cplx(R1, I1));
//				cout2 << "psi(site " << site_number << ", up) = " <<
//				 	eltC(psi(site_number), spin_ind = 1, ri = 1) << "\n";
//				cout2 << "psi(site " << site_number << ", down) = " <<
//				 	eltC(psi(site_number), spin_ind = 2, ri = 1) << "\n";
			}
			else if (site_number == N)
			{
				Index li = leftLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, li = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, li = 1, Cplx(R1, I1));
//				cout2 << "psi(site " << site_number << ", up) = " <<
//				 	eltC(psi(site_number), spin_ind = 1, li = 1) << "\n";
//				cout2 << "psi(site " << site_number << ", down) = " <<
//				 	eltC(psi(site_number), spin_ind = 2, li = 1) << "\n";
			}
			else
			{
				Index li = leftLinkIndex(psi, site_number);
				Index ri = rightLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, li = 1, ri = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, li = 1, ri = 1, Cplx(R1, I1));
//				cout2 << "psi(site " << site_number << ", up) = " <<
//				 	eltC(psi(site_number), spin_ind = 1, li = 1, ri = 1) << "\n";
//				cout2 << "psi(site " << site_number << ", down) = " <<
//				 	eltC(psi(site_number), spin_ind = 2, li = 1, ri = 1) << "\n";
			}
		}
		psi.orthogonalize(Args("Cutoff", 1e-6));

		psi_defined = true;
		//Compute the density matrix rho associated to the pure state |psi>
		cout2 << "Constructing rho from |psi>... Norm of |psi>: (<psi|psi>)^(1/2)=" << norm(psi) <<" ";
  		cout2.flush();
		C.psi2rho(psi, argsRho);
		cout2 << "psi2rho done.\nMax bond dimension of rho:" << maxLinkDim(C.rho)<<"\n";
		cout2.flush();
		for (int site_number = 1; site_number <= N; site_number++)
        {
        	const unsigned int i_site = (unsigned int)(site_number - 1);
        	double b = a_mixed_state[i_site];
        	if (b >= 0. && b <= 1. )
        	{
				double d_i=di_mixed_state[i_site];
				double d_r=dr_mixed_state[i_site];
				psi_defined=false;
				auto pauli_ind = siteIndex(C.rho, site_number);
				if (site_number == 1) {
						  Index ri = rightLinkIndex(C.rho, site_number);
						  int dr=ri.dim();
						  for (int r=1;r<=dr;r++) {
								complex<double> trace=0;
								for (int p=1;p<=4;p+=3)
    								trace+=eltC(C.rho.ref(site_number), ri = r, pauli_ind = p);
								C.rho.ref(site_number).set(pauli_ind = 1, ri = r, b * trace); // |u><u|
						  		C.rho.ref(site_number).set(pauli_ind = 2, ri = r, complex<double>(d_r,-d_i)* trace); // |d><u|
						  		C.rho.ref(site_number).set(pauli_ind = 3, ri = r, complex<double>(d_r,d_i)* trace); // |u><d|
						  		C.rho.ref(site_number).set(pauli_ind = 4, ri = r, (1. - b) * trace); // |d><d|
				 		   }
				}
				else if (site_number == N) {
  					  Index li = leftLinkIndex(C.rho, site_number);
					  int dl=li.dim();
					  for (int l=1;l<=dl;l++) {
							complex<double> trace=0;
							for (int p=1;p<=4;p+=3) {
								trace+=eltC(C.rho.ref(site_number),li=l,pauli_ind=p);
							}
							C.rho.ref(site_number).set(pauli_ind = 1, li = l, b*trace);// |u><u|
							C.rho.ref(site_number).set(pauli_ind = 2, li = l, complex<double>(d_r,-d_i)* trace);// |d><u|
							C.rho.ref(site_number).set(pauli_ind = 3, li = l, complex<double>(d_r,d_i)* trace);// |u><d|
							C.rho.ref(site_number).set(pauli_ind = 4, li = l, (1. - b)*trace);// |d><d|
					   }
				}
				else {
						  Index li = leftLinkIndex(C.rho, site_number);
						  Index ri = rightLinkIndex(C.rho, site_number);
						  int dl=li.dim();int dr=ri.dim();
						   for (int r=1;r<=dr;r++)
							for (int l=1;l<=dl;l++) {
								complex<double> trace=0;
								for (int p=1;p<=4;p+=3) {
    								trace+=eltC(C.rho.ref(site_number),li=l,ri=r,pauli_ind=p);
								}
								C.rho.ref(site_number).set(pauli_ind = 1, li = l, ri = r, b*trace);// |u><u|
						  		C.rho.ref(site_number).set(pauli_ind = 2, li = l, ri = r, complex<double>(d_r,-d_i)* trace);// |d><u|
						  		C.rho.ref(site_number).set(pauli_ind = 3, li = l, ri = r, complex<double>(d_r,d_i)* trace);// |u><d|
						  		C.rho.ref(site_number).set(pauli_ind = 4, li = l, ri = r, (1. - b)*trace);// |d><d|
				 			}
				}
        	}
		}

	}

	if (b_cz_pairs) {
			psi_defined=false;
			cout2 << "Application of CZ gates on the requested " << init_cz_gates.size() / 2 << " pairs.\n"; cout2.flush();
			for (unsigned int n = 0; n < init_cz_gates.size(); n += 2)
			{
				const int i = init_cz_gates[n], j = init_cz_gates[n + 1];
				const unsigned int i0 = (unsigned int)(i - 1);
				const unsigned int j0 = (unsigned int)(j - 1);
				if (param.boolval("b_apply_gate_compression"))
						ApplyControlledZGate(C.rho,C.siteops,i,j,argsRho);
					else
						ApplyControlledZGate(C.rho,C.siteops,i,j);
			}
	}
	// If required, finish the initial state construction by applying some single-qubit gates
	vector<double> gate_times;
	vector<string> gate_names;
	vector<unsigned int> gate_i;
	vector<unsigned int> gate_j;
	vector<string> apply_gates = param.stringvec("apply_gates");
	for (unsigned int n=0;n<apply_gates.size();n++) {
		vector<string> vs=split(apply_gates[n],' ');
		if (vs.size()!=3 && vs.size()!=4)
			cout2 << "Error: expecting time gate qubit1 (qubit2) but got '" << apply_gates[n] << "'.\n", exit(1);
		double time=0;
		int i=0,j=0;
		try {
			time=stod(vs[0]);
		}
		catch (...) {
			cout2 << "Error: " << vs[0] << " is not a double (expecting a time value) in '"<<apply_gates[n] << "'.\n", exit(1);
		}
		if ((time - t_0) < (-.1 * tau) || (time - t_f) > .1 * tau)
			cout2 << "Error: time " << vs[0] << " defined in parameter apply_gates is not between t_init and t_final.\n", exit(1);
//		if (fmod(abs(time / tau), 1.) > 0.1 && fmod(abs(time / tau), 1.) < 0.9)
//			cout2 << "Error: time " << vs[0] << " defined in parameter apply_gates is not close to an integer multiple of tau.\n", exit(1);
		if (fmod(abs(time / tau), 1.) == 0.5)
			time += 0.01 * tau; // if requested time is exactly at 0.5 between two tau steps,
			// shift it a bit up in order for the simple rounding below to result in a unique gate application

		string sgate = vs[1];
		transform(sgate.begin(), sgate.end(), sgate.begin(), ::toupper);
		if (sgate=="X" || sgate=="Y" || sgate=="Z" || sgate=="SX" || sgate=="H" || sgate=="U" || sgate=="D") {
				try {
					i=stod(vs[2]);
					if (i<1 || i>N) cout2 << "Error: qubit index "<<i<<" out of range in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				catch (...) {
					cout2 << "Error: " << vs[2] << " is not a valid integer (expecting a qubit number) in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				if (vs.size()==4)
					cout2 << "Error: too many arguments for gate " << sgate <<" in '"<<apply_gates[n] << "'.\n", exit(1);
				gate_times.push_back(time);
				gate_names.push_back(sgate);
				gate_i.push_back(i);
				gate_j.push_back(-1);
		} else if (sgate=="CZ" || sgate=="CX") {
				if (vs.size()==3)
					cout2 << "Error: missing argument for gate " << sgate <<" in '"<<apply_gates[n] << "'.\n", exit(1);
				try {
					i=stod(vs[2]);
					if (i<1 || i>N) cout2 << "Error: qubit index "<<i<<" out of range in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				catch (...) {
					cout2 << "Error: " << vs[2] << " is not a valid integer (expecting a qubit number) in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				try {
					j=stod(vs[3]);
					if (j<1 || j>N) cout2 << "Error: qubit index "<<j<<" out of range in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				catch (...) {
					cout2 << "Error: " << vs[3] << " is not a valid integer (expecting a qubit number) in '"<<apply_gates[n] << "'.\n", exit(1);
				}
				gate_times.push_back(time);
				gate_names.push_back(sgate);
				gate_i.push_back(i);
				gate_j.push_back(j);
		} else {
			cout2 << "Error: "<<sgate<<" is not a valid gate name.\n", exit(1);
		}
	}

	vector<string> custom_obs = param.stringvec("custom_observables", ';');
	const bool b_custom_obs = custom_obs.size() > 0;
	// each of custom_obs's strings now corresponds to one observable. The string can have two formats:
	// "obs_name obs_type:gate_name q0 q1, ..." with obs_type == 'g' or:
	// "obs_name obs_type:op_name q0, ..." with obs_type == 'o'.
	// Therefore it should first be split using the ':' delimiter (to get two sub-strings), and then
	// the first resulting string split by ' ', and the second resulting string should be split by ','.
    // obs_name is the name to write in the output file
    // obs_type is 'g' for gates and 'o' for operators
    // gate_name is similar to the names in 'apply_gates', and q0 q1 are the qubits (q1 for 2Q gates only).

	vector<MPS> ProjectorList;
	vector<string> ProjectorNames;
	vector<string> OperatorObsNames;
	vector<vector<string>> OperatorObs;
	vector<vector<int>> OperatorObsQubits;
	for (string c_obs : custom_obs) {
		vector<string> obs_defs=split(c_obs, ':');
		if (obs_defs.size() != 2)
		    cout2 << "No gate or operator data in parameter custom_observables: " << c_obs << "\n", exit(0);
    	vector<string> obs_head=split(obs_defs[0], ' ');
        if (obs_head[1] == "g")
        {
            // cout2 << "Custom observable defined by gates, name: " << obs_head[0] << ", gates: "<< obs_defs[1] << ";  ";
            ProjectorNames.push_back(obs_head[0]);
            auto psi0_ini = InitState(C.sites);
            for (int i = 1; i <= N; ++i) psi0_ini.set(i, "Up");// Start with all spins up
            MPS psi0 = MPS(psi0_ini);
            ApplyListOfGatesOnAPureState(obs_defs[1],psi0,C);
            psi0.position(1);
            ProjectorList.push_back(MPS(C.siteops));
            MPS& proj=ProjectorList.back();
            C.psi2rho(psi0,proj);
            int max_bd = maxLinkDim(proj);
            cout2 << "MaxBondDim(gates|0...0><0...0|gates) = " << max_bd << "\n";
        }
        else if (obs_head[1] == "o")
        {
            // cout2 << "Custom observable defined by operators, name: " << obs_head[0] << ", operators: "<< obs_defs[1] << ";  ";
            OperatorObsNames.push_back(obs_head[0]);
            vector<string> obs_ops = vector<string>();
            vector<int> obs_qubits = vector<int>();
            StringToOperatorsList(obs_defs[1], obs_ops, obs_qubits);
            OperatorObs.push_back(obs_ops);
            OperatorObsQubits.push_back(obs_qubits);
        }
        else
		    cout2 << "Type of an observable in custom_observables is unknown (must be 'g' or 'o'): " << c_obs << "\n", exit(0);
	}

	vector<string> collapse = param.stringvec("collapse", ';');
	vector<string> CollapseOpsNames;
	vector<vector<string>> CollapseOps;
	vector<vector<int>> CollapseOpsQubits;
	for (string coll : collapse) {
		vector<string> op_defs=split(coll, ':');
		if (op_defs.size() != 2)
		    cout2 << "No operator data in parameter collapse: " << coll << "\n", exit(0);
    	vector<string> op_head=split(op_defs[0], ' ');
        if (op_head[1] == "o")
        {
            // cout2 << "Collapse operator, name: " << op_head[0] << ", operators: "<< op_defs[1] << ";  ";
            CollapseOpsNames.push_back(op_head[0]);
            vector<string> coll_ops = vector<string>();
            vector<int> coll_qubits = vector<int>();
            StringToOperatorsList(op_defs[1], coll_ops, coll_qubits);
            CollapseOps.push_back(coll_ops);
            CollapseOpsQubits.push_back(coll_qubits);
        }
        else
		    cout2 << "Type of an operator in collapse is unknown (must be 'o'): " << coll << "\n", exit(0);
	}
    cout2.flush();

	if (param.val("b_initial_rho_compression") != 0) {

		cout2 << "C.rho.orthogonalize...";
		cout2.flush();
		C.rho.orthogonalize(Args("Cutoff", param.val("cut_off_rho"), "MaxDim", param.longval("max_dim_rho")));
		cout2 << "done.\n";
		cout2 << "New max bond dimension of rho:" << maxLinkDim(C.rho) << "\n";
		cout2.flush();
	}

	Cplx tr = C.trace_rho();
	cout2 << "Tr{rho} before re-normalization: " << tr << "\n"; cout2.flush();
	C.rho /= tr; //Normalize rho so that Tr[rho]=1
	tr = C.trace_rho();
	cout2 << "Tr{rho} after re-normalization: " << tr << "\n"; cout2.flush();

	cout2 << "Tr{rho^2} ="; cout2.flush();
	const Cplx tr2 = C.trace_rho2();
	cout2 << tr2 << "\n"; cout2.flush();

	if (psi_defined)
	{
		if (std::abs(tr - 1) > TRACE_RHO_THRESHOLD || std::abs(tr2 - 1) > TRACE_RHO_THRESHOLD)
			cout2 << "Error, those traces should be 1 for a pure state |psi><psi|.\n", C.rho /= tr;
		//Check a few simple observables, using rho and psi
		vector<string> ops = {"Sz", "S+", "S-", "Sx", "Sy"};
		vector<double> psi_factor = {2., 1., 1., 2., 2.};
		double n2=norm(psi);n2*=n2;//n2=<psi|psi>
		double err = 0;
		for (int i = 1; i <= N; i++)
		{
			for (unsigned int o = 0; o < ops.size(); o++)
			{
				const string &opname = ops[o];
				Cplx with_rho = C.Expect(opname, i);
				Cplx with_psi = psi_factor[o] * PureStateObs(opname, psi, i, C.sites)/n2;
				err += std::abs(with_rho - with_psi);
				if (std::abs(with_rho - with_psi) > TRACE_RHO_THRESHOLD)
				cout2 << "Error: <psi|" << opname << "(" << i << ")|psi> / <psi|psi> =" << with_psi << "\t"
					<< "Tr[rho*" << opname << "(" << i << ")]=" << with_rho << "\n",
					exit(1);
			}
		}
	}

	//-----------------------------------------------------
	//Construct the Lindbladian from the parameters (unitary and dissipative terms)

	SetLindbladian(C, param, lattice);

	//-----------------------------------------------------
	cout2 << "Compute exp(i*tau*L) as an MPO... ";
	cout2.flush();

	const int o = param.val("trotter_order");
	TimeEvolver TE; //Object defined in "TimeEvolution.h" and "TimeEvolution.cc"
	TE.init(tau, C.Lindbladian, argsRho, o);
	cout2 << "done.\n";
	cout2.flush();
	cout2 << "Largest bond dimension of the MPO exp(tau*L): " << maxLinkDim(TE.expL1) << "\n";
	const int n_steps = int(t_total / tau);

	// Open output files
	ofstream file_1q, file_2q, file_3q, file_global, file_custom, file_collapse;
	file_global.open(output_prefix + ".global.dat"); // Always written to.
	file_global.precision(15);
	file_global << "#time\tquantity\tvalue" << endl;
	if (b_custom_obs)
	{
        file_custom.open(output_prefix + ".obs-cu.dat");
        file_custom.precision(15);
        file_custom << "#time\tobservable\tvalue" << endl;
    }
    if (collapse.size() && !b_custom_obs)
    {
        cout2 << "Error: collapse operators list is nonempty, but no custom observables defined.\n", exit(1);
    }
	// Some preparation/checks for the 1-qubit observables
	auto components = param.stringvec("1q_components");
	vector<long> sit = param.longvec("1q_indices");
	if (components.size())
	{
        for (auto &s : components)
        {
            if (s.length() != 1)
                cout2 << "Error: " << s << " is an unknown 1-qubit component (should be in {x,y,z} or in {X,Y,Z}).\n", exit(1);
            char c = toupper(s[0]);
            if (c != 'X' && c != 'Y' && c != 'Z')
                cout2 << "Error: " << s << " is an unknown 1-qubit component (should be in {x,y,z} or in {X,Y,Z}).\n", exit(1);
        }
        if (sit.size() == 0)
        { //If no sites are given explicitly we consider all: 1,...,N
            sit.resize(N);
            iota(sit.begin(), sit.end(), 1);
        }
        for (int i : sit)
        {
            if (i < 1 || i > N)
                cout2 << "Error: invalid index i=" << i << " found in list `1q_indices`.\n", exit(1);
        }
        file_1q.open(output_prefix + ".obs-1q.dat");
        file_1q.precision(15);
        file_1q << "#time\toperator\tindex\tvalue" << endl;
    }

	// Some preparation/checks for the 2-qubit observables
	auto components2 = param.stringvec("2q_components");
	vector<long> sit2 = param.longvec("2q_indices");
	if (components2.size())
	{
        if (sit2.size() % 2 == 1)
            cout2 << "Error: the list of indices given in the parameter `2q_indices` should have an even length.\n", exit(1);
        if (sit2.size() == 0)
        { //If no sites are given explicitly we consider all pairs 1,2,1,3,...,1,N,    2,1,2,3,2,4,...,2,N,  ...  N,N-1
            for (int i = 1; i <= N; i++)
                for (int j = 1; j <= N; j++)
                {
                    if (i != j)
                        sit2.push_back(i), sit2.push_back(j);
                }
        }
        validate_2q_list(sit2, N, "q2_indices");
        for (auto &s : components2)
        {
            if (s.length() != 2)
                cout2 << "Error: " << s << " is an unknown 2-qubit component (should be a pair in (x,y,z)*(x,y,z)).\n", exit(1);
            for (int n = 0; n <= 1; n++)
            {
                char c = toupper(s[n]);
                if (c != 'X' && c != 'Y' && c != 'Z')
                    cout2 << "Error: " << s << " is an unknown component (should be a pair in (x,y,z)*(x,y,z)).\n", exit(1);
            }
        }
        file_2q.open(output_prefix + ".obs-2q.dat");
        file_2q.precision(15);
        file_2q << "#time\toperator\tindex_1\tindex_2\tvalue" << endl;
	}

	// Some preparation/checks for the 3-qubit observables
	auto components3 = param.stringvec("3q_components");
	vector<long> sit3 = param.longvec("3q_indices");
	if (components3.size())
	{
        if (sit3.size() % 3 > 0)
            cout2 << "Error: the list of indices given in the parameter `3q_indices` should be multiple of three.\n", exit(1);
        if (sit3.size() == 0)
        { //If no sites are given explicitly we consider all triples
            for (int i = 1; i <= N; i++)
                for (int j = 1; j <= N; j++)
                    for (int k = 1; k <= N; k++)
                    {
                        if (i != j && i != k && j != k)
                            sit3.push_back(i), sit3.push_back(j), sit3.push_back(k);
                    }
        }
        validate_3q_list(sit3, N, "q3_indices");
        for (auto &s : components3)
        {
            if (s.length() != 3)
                cout2 << "Error: " << s << " is an unknown 3-qubit component (should be a triplet in (x,y,z)*(x,y,z)*(x,y,z)).\n", exit(1);
            for (int n = 0; n <= 2; n++)
            {
                char c = toupper(s[n]);
                if (c != 'X' && c != 'Y' && c != 'Z')
                    cout2 << "Error: " << s << " is an unknown component (should be a triplet in (x,y,z)*(x,y,z)*(x,y,z)).\n", exit(1);
            }
        }
        file_3q.open(output_prefix + ".obs-3q.dat");
        file_3q.precision(15);
        file_3q << "#time\toperator\tindex_1\tindex_2\tindex_3\tvalue" << endl;
    }

	//-----------------------------------------------------
	const int output_step = param.longval("output_step");
	auto t_init_end = steady_clock::now();
	auto duration_ms = duration_cast<milliseconds>(t_init_end - t_start_sim);
	cout2 << "\nSimulation initialization took: " << duration_ms.count() / 1000. << "s" << "\n";
	cout2.flush();

	char buf[100];
	const long force_rho_Hermitian_step = param.longval("force_rho_hermitian_step");
 	const bool b_quiet = param.boolval("b_quiet");
 	cout2.quiet(b_quiet);
 	double t = t_0;
	for (int n = 0; n <= n_steps; n++)
	{
		t = t_0 + n * tau;
		// Print data about this time step
		auto t_now = steady_clock::now();
		auto tot_duration = duration_cast<milliseconds>(t_now - t_start_sim);
		sprintf(buf, "%.2fhr", tot_duration.count() / 3600000.);
		cout2 << "\nSolution time t = " << t << " ----------------------";
		cout2 << " Total run duration: " << buf << "\n";
		cout2.flush();

		//If the time corresponds a step where some gates should be applied:
		//note: if several gates are associated to the same time, the application will follow
		// the order of the arguments of 'apply_gate'
		bool b_normalize = false; // Will be set to true if projectors applied to rho.
		for (unsigned int k=0;k<gate_times.size();k++) {
			if (abs(gate_times[k] - t) < (tau / 2.)) {
				// This will make the gate execute at t nearest to the requested time.
				// Earlier there is a validation making sure this is unique.
				if (gate_names[k]=="CZ") {
					cout2<<"\tApplication of gate "<<gate_names[k]<<"("<<gate_i[k]<<","<<gate_j[k]<<")\n";
					if (param.boolval("b_apply_gate_compression"))
						ApplyControlledZGate(C.rho,C.siteops,gate_i[k],gate_j[k],argsRho);
					else
						ApplyControlledZGate(C.rho,C.siteops,gate_i[k],gate_j[k]);
				}
				else if (gate_names[k]=="CX") {
					cout2<<"\tApplication of gate "<<gate_names[k]<<"("<<gate_i[k]<<","<<gate_j[k]<<")\n";
					if (param.boolval("b_apply_gate_compression"))
						ApplyCNOTGate(C.rho,C.siteops,gate_i[k],gate_j[k],argsRho);
					else
						ApplyCNOTGate(C.rho,C.siteops,gate_i[k],gate_j[k]);
				} else {
					if (gate_names[k]=="U") {
      					cout2<<"\tApplication of projector Up" <<"("<<gate_i[k]<<")\n";
					    ApplyProjUp(C.rho,C.siteops,gate_i[k]);
					    b_normalize = true;
					}
					else {
                        if (gate_names[k]=="D") {
          					cout2<<"\tApplication of projector Down" <<"("<<gate_i[k]<<")\n";
                            ApplyProjDn(C.rho,C.siteops,gate_i[k]);
                            b_normalize = true;
                        }
                        else {
        					cout2<<"\tApplication of gate "<<gate_names[k]<<"("<<gate_i[k]<<")\n";
                            if (gate_names[k]=="X") ApplyXGate(C.rho,C.siteops,gate_i[k]);
                            if (gate_names[k]=="Y") ApplyYGate(C.rho,C.siteops,gate_i[k]);
                            if (gate_names[k]=="Z") ApplyZGate(C.rho,C.siteops,gate_i[k]);
                            if (gate_names[k]=="SX") ApplySqrtXGate(C.rho,C.siteops,gate_i[k]);
                            if (gate_names[k]=="H") ApplyHGate(C.rho,C.siteops,gate_i[k]);
                        }
                    }
				}
			}
		}
        if (b_normalize) {
			Cplx z = C.trace_rho();
			cout2 << "\tNormalizing rho after projections applied. Tr{rho}: " << z << "\n";
			if (std::abs(z) < _2_N)
				cout2 << "\t\tNote: this is smaller than 2^(-N)!" << "\n";
				// TODO: This is a somewhat arbitrary threshold for the warning.
    		C.rho /= z;
        }
		cout2.flush();

		if (force_rho_Hermitian_step && (n % force_rho_Hermitian_step) == 0)
			C.MakeRhoHermitian(argsRho);
		if (output_step > 0)
		{
			if ((n % output_step) == 0 || n == n_steps)
			{
				// Print and save output data at initial time, final time, and every output_step time steps

				tr = C.trace_rho();
				const Cplx tr2 = C.trace_rho2();
				const Real osee = OSEE(C.rho, N / 2);
				// This is the operator-space entanglement entropy (OSEE) of rho, associated
				// to a cut in the middle of the system (at the center bond)
				const double S_2 = 1.0 / (1.0 - 2.0) * log(tr2.real());
				const int bd = BondDim(C.rho, N / 2), bd_max = maxLinkDim(C.rho);

				cout2 << "\tTr{rho}: " << tr << ", Renyi Entropy S_2: " << S_2; // << ",\tTr{rho^2} =" << tr2
//				     << "\n\tCenter bond dimension: " << bd << ", Max bond dimension: " << bd_max
				if (!force_rho_Hermitian_step || (n % force_rho_Hermitian_step) != 0)
				     cout2 << "\n\tMax bond dimension: " << bd_max;
				cout2 << "\n\tOperator space entanglement entropy at center bond: " << osee;

//				file_global << t << " \t" << tr.real() << " \t" << S_2 << " \t" << osee << " \t" << bd_max << " \t" << tot_duration.count() << endl;
				file_global << t << " \t" << "tr_rho\t" << tr.real() << "\n";
				file_global << t << " \t" << "S_2\t" << S_2 << "\n";
				file_global << t << " \t" << "OSEE_center\t" << osee << "\n";
				file_global << t << " \t" << "max_bond_dim\t" << bd_max << "\n";
				file_global << t << " \t" << "duration_ms\t" << tot_duration.count() << "\n";
				file_global << endl; // Skip a line between time steps

				// --------------------------------------------------
				// Compute 1-qubit observables and write them to file
				int count = 0;
				auto t_1q_start = steady_clock::now();
                if (components.size())
                {
                    for (long &i : sit)
                    {
                        for (auto &s : components)
                        {
                            string c1("S");
                            c1 += char(tolower(s[0]));
                            Cplx expectation_value = C.Expect(c1, i);
                            if (abs(expectation_value.imag()) > IMAGINARY_THRESHOLD)
                                cout2 << "\tWarning: <S^" << s << "(" << i << ")> = " << expectation_value <<
                                    "; it should be real, but has an imaginary part > " <<
                                    IMAGINARY_THRESHOLD << ".\n";
                            file_1q << t << "\t" << char(toupper(s[0])) << "\t" << i
                                  << "\t" << expectation_value.real() << endl;
                            count++;
                        }
    //					file_1q << endl;
                    }
                    file_1q << endl; // Skip a line between time steps
                }
                auto t_1q_end = steady_clock::now();
                if (count)
                {
                    duration_ms = duration_cast<milliseconds>(t_1q_end - t_1q_start);
                    cout2 << "\n\t" << count << " 1-qubit expectation values saved to file. Duration: " << duration_ms.count() / 1000. << "s";
    				count = 0;
                }
		        // --------------------------------------------------
				// Compute 2-qubit observables and write them to file
				if (components2.size())
				{
                    for (unsigned int n = 0; n < sit2.size(); n += 2)
                    {
                        const int i = sit2[n], j = sit2[n + 1];
                        //Loop over components
                        for (auto &s : components2)
                        {
                            string c1("S"), c2("S");
                            c1 += char(tolower(s[0]));
                            c2 += char(tolower(s[1]));
                            Cplx expectation_value = C.Expect(c1, i, c2, j);
                            if (abs(expectation_value.imag()) > IMAGINARY_THRESHOLD)
                                cout2 << "\tWarning: <" << c1 << "(" << i << ")" << c2 << "(" << j <<
                                    ")> = " << expectation_value <<
                                    "; it should be real, but has an imaginary part > " <<
                                    IMAGINARY_THRESHOLD << ".\n";
                            file_2q << t << "\t" << char(toupper(s[0])) << char(toupper(s[1])) << "\t" << i << "\t" << j << "\t" << expectation_value.real() << endl;
                            count++;
                        }
                    }
                    file_2q << endl; //Skip a line between time steps
                }
				auto t_2q_end = steady_clock::now();
				if (count)
				{
					duration_ms = duration_cast<milliseconds>(t_2q_end - t_1q_end);
					cout2 << "\n\t" << count << " 2-qubit expectation values saved to file. Duration: " << duration_ms.count() / 1000. << "s";
    				count = 0;
				}
				// --------------------------------------------------
				// Compute 3-qubit observables and write them to file
				if (components3.size())
                {
                    for (unsigned int n = 0; n < sit3.size(); n += 3)
                    {
                        const int i = sit3[n], j = sit3[n + 1], k = sit3[n + 2];
                        //Loop over components
                        for (auto &s : components3)
                        {
                            string c1("S"), c2("S"),c3("S");
                            c1 += char(tolower(s[0]));
                            c2 += char(tolower(s[1]));
                            c3 += char(tolower(s[2]));
                            Cplx expectation_value = C.Expect(c1, i, c2, j,c3,k);
                            if (abs(expectation_value.imag()) > IMAGINARY_THRESHOLD)
                                cout2 << "\tWarning: <" << c1 << "(" << i << ")" << c2 << "(" << j <<
                                    ")" << c3 << "(" << k <<")" << expectation_value <<
                                    "; it should be real, but has an imaginary part > " <<
                                    IMAGINARY_THRESHOLD << ".\n";
                            file_3q << t << "\t" << char(toupper(s[0])) << char(toupper(s[1])) << char(toupper(s[2]))
                            << "\t" << i << "\t" << j << "\t" << k <<"\t" << expectation_value.real() << endl;
                            count++;
                        }
                    }
                    file_3q << endl; //Skip a line between time steps
                }
				auto t_3q_end = steady_clock::now();
				if (count)
                {
                    duration_ms = duration_cast<milliseconds>(t_3q_end - t_2q_end);
                    cout2 << "\n\t" << count << " 3-qubit expectation values saved to file. Duration: " << duration_ms.count() / 1000. << "s";
    				count = 0;
                }

				// --------------------------------------------------
				// Custom observables
				if (b_custom_obs)
				{
                    if (ProjectorList.size()) {
                        int c=0;
                        for (MPS& proj:ProjectorList) {
                            Cplx op_val = innerC(proj,C.rho);
                            file_custom << t << " \t" << ProjectorNames[c] << "\t" << op_val.real() << endl;
                            c++;
                        }
                        count += c;
                    }
                    if (OperatorObsNames.size()) {
                        int c=0;
                        for (string& s_op_obs : OperatorObsNames)
                        {
                            Cplx op_val = C.Expect(OperatorObs[c], OperatorObsQubits[c]);
                            file_custom << t << " \t" << s_op_obs << "\t" << op_val.real() << endl;
                            c++;
                        }
                        count += c;
                    }
                    file_custom << endl; //Skip a line between time steps
                }
				auto t_cu_end = steady_clock::now();
				if (count)
                {
                    duration_ms = duration_cast<milliseconds>(t_cu_end - t_3q_end);
                    cout2 << "\n\t" << count << " custom expectation values saved to file. Duration: " << duration_ms.count() / 1000. << "s";
                    count = 0;
                }
				cout2 << "\n";
				cout2.flush();
			}
		}
		if (n < n_steps)
		{
			cout2 << "\t" << "Time evolving the state -> ";
			cout2.flush();
			auto t_evolve_start = steady_clock::now();
			TE.evolve(C.rho);
			auto t_evolve_end = steady_clock::now();
			duration_ms = duration_cast<milliseconds>(t_evolve_end - t_evolve_start);
			cout2 << "done. Duration: " << duration_ms.count() / 1000. << "s" << "\n";

			Cplx z = C.trace_rho();
			if (std::abs(z - 1) > TRACE_RHO_THRESHOLD)
				cout2 << "\nWarning: Tr[rho] != 1 :" << z << "\n";
			if (param.val("b_force_rho_trace") != 0)
				C.rho /= z;
			cout2.flush();
		}
	}
	cout2.quiet(false);
	cout2 << "\nSimulation ended.\n";

	bool b_save_state = param.boolval("b_save_final_state");
	if (b_save_state)
	{
		string f1 = output_prefix + ".state.ops";
		writeToFile(f1, C.siteops);
		string f2 = output_prefix + ".state.rho";
		writeToFile(f2, C.rho);
		string f3 = output_prefix + ".state.sites";
		writeToFile(f3, C.sites);
		cout2 << "The final state was saved to disk, using 3 files:\n" << f1 << "\n" << f2 << "\n" << f3 << "\n";
	}
	if (file_1q.is_open())
    	file_1q.close();
	if (file_2q.is_open())
        file_2q.close();
	if (file_3q.is_open())
    	file_3q.close();
	if (file_global.is_open())
    	file_global.close();
	if (file_custom.is_open())
    	file_custom.close();

    if (collapse.size())
    {
        // auto t_collapse_start = steady_clock::now();
    	cout2 << "\nStarting evaluation of " << collapse.size() << " collapse operators.\n";
        file_collapse.open(output_prefix + ".obs-co.dat");
        file_collapse.precision(15);
        file_collapse << "#time\tcollapse_op\tobservable\tvalue" << endl;

        int i_coll = 0;
        MPS rho_0(C.rho);  // Keep a copy of the uncollapsed state
        for (string& s_coll_op : CollapseOpsNames)
        {
            // collapse rho and send it to the observables
            auto t_cu_start = steady_clock::now();
            vector<string> &op_names = CollapseOps[i_coll];
            vector<int> &op_qubits = CollapseOpsQubits[i_coll];
            int i_op = 0 ;
            for (string &op_name: op_names)
            {
                if (op_name[1] == 'u') {
                    ApplyProjUp(C.rho,C.siteops,op_qubits[i_op]);
                }
                else if (op_name[1]=='d') {
                    ApplyProjDn(C.rho,C.siteops,op_qubits[i_op]);
                }
                else
                    cout2 << "Error: Collapse operator " << s_coll_op <<
                        " can only contain standard-basis projectors (u or d)," <<
                        " but contains " << op_name << ".\n", exit(1);
                i_op++;
            }
            Cplx z = C.trace_rho();
			cout2 << "\tNormalizing rho after collapse projector applied. Tr{rho}: " << z << "\n";
			if (std::abs(z) < _2_N)
				cout2 << "\t\tNote: this is smaller than 2^(-N)!" << "\n";
				// TODO: This is a somewhat arbitrary threshold for the warning.
       		C.rho /= z;

            int count = 0;
            if (ProjectorList.size()) {
                int c=0;
                for (MPS& proj:ProjectorList) {
                    Cplx op_val = innerC(proj,C.rho);
                    file_collapse << t << " \t" << s_coll_op << " \t" << ProjectorNames[c] << "\t" << op_val.real() << endl;
                    c++;
                }
                count += c;
            }
            if (OperatorObsNames.size()) {
                int c=0;
                for (string& s_op_obs : OperatorObsNames)
                {
                    Cplx op_val = C.Expect(OperatorObs[c], OperatorObsQubits[c]);
                    file_collapse << t << " \t" << s_coll_op << " \t" << s_op_obs << "\t" << op_val.real() << endl;
                    c++;
                }
                count += c;
            }
            file_collapse << endl; //Skip a line between time steps
            C.rho = MPS(rho_0);
            auto t_cu_end = steady_clock::now();
            if (count)
            {
                duration_ms = duration_cast<milliseconds>(t_cu_end - t_cu_start);
                cout2 << "\t" << count << " custom expectation values saved to file. Duration: " << duration_ms.count() / 1000. << "s\n";
                count = 0;
            }
            i_coll++;
        }

    	file_collapse.close();
        // auto t_collapse_end = steady_clock::now();
    }

	auto t_end_sim = steady_clock::now();
	auto tot_duration = duration_cast<seconds>(t_end_sim - t_start_sim);
	sprintf(buf, "%.2fhr", tot_duration.count() / 3600.);
	cout2 << "\nTotal simulation duration: " << buf << "\n";
	cout2.flush();
	log_file.close();
	return 0;
}

void validate_2q_list(vector<long> &vect, int N, string const &list_name)
{
	for (unsigned int n = 0; n < vect.size(); n += 2)
	{
		const int i = vect[n], j = vect[n + 1];
		if (i < 1 || i > N)
			cout2 << "Error: invalid index i =" << i << " found in list `" << list_name <<
				"`.\n", exit(1);
		if (j < 1 || j > N)
			cout2 << "Error: invalid index i =" << j << " found in list `" << list_name <<
				"`.\n", exit(1);
		if (i == j)
			cout2 << "Error: an invalid identical index pair (" << i << ") found in list `" <<
				list_name << "`. Two-qubit operators must involve two distinct qubits.\n", exit(1);
	}
}

void validate_3q_list(vector<long> &vect, int N, string const &list_name)
{
	for (unsigned int n = 0; n < vect.size(); n += 3)
	{
		const int i = vect[n], j = vect[n + 1], k=vect[n+2];
		if (i < 1 || i > N)
			cout2 << "Error: invalid index i =" << i << " found in list `" << list_name <<
				"`.\n", exit(1);
		if (j < 1 || j > N)
			cout2 << "Error: invalid index j =" << j << " found in list `" << list_name <<
				"`.\n", exit(1);
		if (k < 1 || k > N)
			cout2 << "Error: invalid index k =" << k << " found in list `" << list_name <<
				"`.\n", exit(1);
		if (i == j || i == k || k == j)
			cout2 << "Error: invalid index triplet (" << i <<","<< j <<","<< k << ") found in list `" <<
				list_name << "`. 3-qubit operators must involve 3 distinct qubits.\n", exit(1);
	}
}


// Old initialization code
/*
  if (param.stringval("load_purestate_file") != "" && param.stringval("load_state_file") != "")
    cout2 << "Error, conflict in parameters:load_purestate_file=" << param.stringval("load_purestate_file")
         << " and load_state_file=" << param.stringval("load_purestate_file") << ". They should not be both defined\n",
        exit(1);

  //-----------------------------------------------------
  //Hamiltonian used to define the initial (pure) state at t=0 (its ground-state)
  // in the examples below (xm_init, ...) H0 is simply a external magnetic field,
  // but more complicated interactions can alo be considered
  AutoMPO ampo(C.sites);

  if (param.val("xm_init") != 0)
  {
    // To initialize the system in a state whera all spins point in the -x direction,
    // we simply construct an Hamiltonian H0 with a magnetic field pointing in the x direction
    for (int j = 1; j <= N; ++j)
      ampo += +1.0, "Sx", j;
  }
  else
  {
    if (param.val("y_init") != 0)
    {
      // To initialize the system in a state whera all spins point in the x direction,
      // we simply construct an Hamiltonian H0 with a magnetic field pointing in the -y direction
      for (int j = 1; j <= N; ++j)
        ampo += -1.0, "Sy", j;
    }
    else
    {
      if (param.val("x_init") != 0)
      {
        // To initialize the system in a state whera all spins point in the x direction,
        // we simply construct an Hamiltonian H0 with a magnetic field pointing in the x direction
        for (int j = 1; j <= N; ++j)
          ampo += -1.0, "Sx", j;
      }
      else
      {
        if (param.val("ym_init") != 0)
        {
          // To initialize the system in a state whera all spins point in the -y direction,
          // we simply construct an Hamiltonian H0 with a magnetic field pointing in the y direction
          for (int j = 1; j <= N; ++j)
            ampo += 1.0, "Sy", j;
        }
      }
    }
  }
  auto H0 = toMPO(ampo);


  if (param.longval("rho_inf_init") == 0)
  {
    if (param.stringval("load_state_file") == "")

    { //Start from a wave-function (pure state)
      if (param.stringval("load_purestate_file") == "")
      {
        // Set the initial wavefunction matrix product state
        auto initState = InitState(C.sites);

        if (param.longval("up_init") == 0 && param.longval("down_init") == 0)
        {
          for (int i = 1; i <= N; ++i)
          { // Start the DMRG from the Neel state (|Neel> = |up|down|up|down|up|...> )
            if (i % 2 == 1)
              initState.set(i, "Dn");
            else
              initState.set(i, "Up");
          }
        }
        else
        {
          if (param.longval("up_init") != 0 && param.longval("down_init") != 0)
            cout2 << "Error: conflicting initialization options for the DMRG:up_init and  down_init.\n", exit(1);

          if (param.longval("up_init") != 0)
          {
            for (int i = 1; i <= N; ++i) // Start the DMRG from all spins up
              initState.set(i, "Up");
          }
          if (param.longval("down_init") != 0)
          {
            for (int i = 1; i <= N; ++i) // Start the DMRG from all spins Down
              initState.set(i, "Dn");
          }
        }
        //Find the ground-state of H0. THis will be the initial state of the time evolution
        psi = MPS(initState);
        auto sweeps = Sweeps(param.longval("sweeps"));
        //Specify max number of states kept each sweep
        const int m = param.longval("max_dim_psi");
        //We specify bvelow how the max. bond dimension should be increased along the DMRG sweeps
        sweeps.maxdim() = min(5, m), min(5, m), min(10, m), min(10, m), min(10, m), min(20, m), min(20, m), min(50, m), min(50, m), min(100, m), min(100, m), min(200, m), min(200, m), m;
        sweeps.cutoff() = param.val("cut_off_psi");
        //Run the DMRG algorithm
        MyDMRGObserver obs(psi, param.val("energy")); //Convergence criterium on the energy passed to the DMRGObserver

        const double energy = dmrg(psi, H0, sweeps, obs, "Quiet");
        cout2 << "Initial energy=" << energy << "\n";
        if (param.stringval("save_purestate_file") != "")
        {
          string file_name = param.stringval("save_purestate_file");
          file_name += "_N=" + to_string(N);
          string f1 = file_name + ".ops";
          writeToFile(f1, C.siteops);
          string f2 = file_name + ".psi";
          writeToFile(f2, C.rho);
          string f3 = file_name + ".sites";
          writeToFile(f3, C.sites);
          writeToFile(f3, C.sites);
          writeToFile(f1, C.siteops);
          writeToFile(f2, psi);
          cout2 << "the final pure state was written to disk, in files " << f1 << ", " << f2 << " and " << f3 << ".\n";
        }
      }
      else
      { //read psi from a file
        string file_name = param.stringval("load_purestate_file");
        file_name += "_N=" + to_string(N) + ".psi";
        psi = MPS(C.sites);
        cout2 << "Read the initial pure state (wave-function) from file '" << file_name << "'...";
        cout2.flush();
        readFromFile(file_name, psi);
        cout2 << "done.\n";
      }
      psi_defined = true;

      //Compute the density matrix rho associated to the pure state |psi>
      C.psi2rho(psi, argsRho);
      cout2 << "psi2rho done.\n";
      cout2.flush();
    }
    else
    { //Read the density matrix from disk

      string file_name = param.stringval("load_state_file");
      file_name += "_N=" + to_string(N) + ".rho";
      cout2 << "Read the initial rho from the file '" << file_name << "'...";
      cout2.flush();
      readFromFile(file_name, C.rho);
      cout2 << "done.\n";
      if (param.val("InitialOrthoRho") != 0)
      {
        cout2 << "C.rho.orthogonalize...";
        cout2.flush();
        C.rho.orthogonalize(Args("Cutoff", param.val("cut_off_rho"), "MaxDim", param.longval("max_dim_rho")));
        cout2 << "done.\n";
        cout2.flush();
      }
    }
  }
  else
  {
    //Start from rho ~ identity (infinite temperature density matrix)
    cout2 << "Initialize rho ~ identity (infinite temperature density matrix).\n";
    C.rho = C.Identity;
  }
*/
