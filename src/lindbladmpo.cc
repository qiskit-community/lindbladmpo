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
#include <chrono>

using namespace itensor;
using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
	ModelParameters param;

	auto start_sim = steady_clock::now();
	// Start tracking simulation time

	// Read the input parameters given in the command line. Default values are substituted for all
	// missing parameters, as defined in the files SimulationParameters.h and ModelParameters.h
	param.ReadArguments(argc, argv);

	// If a filename was given as the parameter `input_file`, read the parameters given in this file.
	string inputfilename = param.stringval("input_file");
	if (inputfilename != "")
		param.ReadFromFile(inputfilename);
	// Now `param` contains all parameters.

	cout.precision(10);
	cout << endl;
	//param.check();
	param.Print(cout);

	int Lx = param.longval("l_x");
	int Ly = param.longval("l_y");
	int N;
	Lattice2d lattice;
	if (Lx < 0 || Ly < 0)
		cerr << "Error, invalid (l_x, l_y)=" << Lx << "," << Ly << ".\n", exit(1);
	if (Lx > 0 && Ly > 0)
	{
		//Square lattice
		N = Lx * Ly;
		lattice = Lattice2d(Lx, Ly, param.boolval("b_periodic_x"), param.boolval("b_periodic_y"));
	}
	else
	{
		//User-defined lattice
		N = param.longval("N");
		if (N < 1)
			cerr << "Error, invalid N=" << N << ".\n", exit(1);
		lattice = Lattice2d(N, param.longvec("A_bond_indices"), param.longvec("B_bond_indices"));
	}

	SpinHalfSystem C(N);

	string load_prefix = param.stringval("load_files_prefix");
	if (load_prefix != "")  // || param.stringval("load_purestate_file") != "")
	{
//      if (param.stringval("load_purestate_file") != "")
//          file_name = param.stringval("load_purestate_file");
		string file_name = load_prefix;
		file_name += ".N=" + to_string(N);
		string f1 = file_name + ".state.ops";
		string f3 = file_name + ".state.sites";
		cout << "Opening '" << f1 << "' and '" << f3 << "'...";
		cout.flush();
		readFromFile(f3, C.sites);
		readFromFile(f1, C.siteops);
		cout << "done.\n";
		C.rho = MPS(C.siteops);
		// Rho is still undefined but its structure is initialized using C.siteops
		C.Lindbladian = AutoMPO(C.siteops);
		C.LindbladianDag = AutoMPO(C.siteops);
	}

    C.ConstructIdentity(); // Construct the density matrix corresponding to infinite temperature

	// Note on the option below: if we  "Normalize=true" (default for fitapplyMPO)
	// we would normalize rho such that Tr[rho^2]=1, (norm of the MPS)
	// which is of course not appropriate for a density matrix.
	auto argsRho = Args();
	argsRho.add("Normalize", false);
	argsRho.add("MinDim", param.longval("min_dim_rho"));
	argsRho.add("MaxDim", param.longval("max_dim_rho"));
	argsRho.add("Cutoff", param.val("cut_off_rho"));

	vector<string> a_init = param.stringvec("init_Pauli_state");
	unsigned int a_init_len = a_init.size();
	if (a_init_len != 1 && int(a_init_len) != N)
		cerr << "Error: the parameter init_Pauli_state has " << a_init_len << " value(s) but 1 or " << N << " value(s) were expected.\n", exit(1);
	if (a_init_len == 1)
		a_init = vector<string>(N, a_init[0]);

	MPS psi;
	bool psi_defined = false;
	if (load_prefix != "")
	{ //Read the density matrix from disk
	    string file_name = load_prefix;
		file_name += ".N=" + to_string(N) + ".state.rho";
		cout << "Read the initial rho from the file '" << file_name << "'...";
		cout.flush();
		readFromFile(file_name, C.rho);
		cout << "done.\n";
		if (param.val("b_initial_rho_orthogonalization") != 0)
		{
			cout << "C.rho.orthogonalize...";
			cout.flush();
			C.rho.orthogonalize(Args("Cutoff", param.val("cut_off_rho"), "MinDim", param.longval("min_dim_rho"),
				"MaxDim", param.longval("max_dim_rho")));
			cout << "done.\n";
			cout.flush();
		}
	}
	else
	{
    // Set the initial wavefunction matrix product state
		auto initState = InitState(C.sites);
		for (int i = 1; i <= N; ++i) // Start with all spins up
		initState.set(i, "Up");
		psi = MPS(initState);
		double sqrt05 = pow(.5, .5);
		for (int site_number = 1; site_number <= N; site_number++)
        {
			string s_init = a_init[site_number - 1];
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

            auto spin_ind = siteIndex(psi, site_number);
            if (site_number == 1)
			{
				Index ri = rightLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, ri = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, ri = 1, Cplx(R1, I1));
				cout << "psi(site " << site_number << ",up)=" << elt(psi(site_number), spin_ind = 1, ri = 1) << endl;
				cout << "psi(site " << site_number << ",down)=" << elt(psi(site_number), spin_ind = 2, ri = 1) << endl;
			}
			if (site_number == N)
			{
				Index li = leftLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, li = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, li = 1, Cplx(R1, I1));
				cout << "psi(site " << site_number << ",up)=" << elt(psi(site_number), spin_ind = 1, li = 1) << endl;
				cout << "psi(site " << site_number << ",down)=" << elt(psi(site_number), spin_ind = 2, li = 1) << endl;
			}
			if (site_number > 1 && site_number < N)
			{
				Index li = leftLinkIndex(psi, site_number);
				Index ri = rightLinkIndex(psi, site_number);
				psi.ref(site_number).set(spin_ind = 1, li = 1, ri = 1, Cplx(R0, I0));
				psi.ref(site_number).set(spin_ind = 2, li = 1, ri = 1, Cplx(R1, I1));
				cout << "psi(site " << site_number << ",up)=" << eltC(psi(site_number), spin_ind = 1, li = 1, ri = 1) << endl;
				cout << "psi(site " << site_number << ",down)=" << eltC(psi(site_number), spin_ind = 2, li = 1, ri = 1) << endl;
			}
		}
		psi.orthogonalize(Args("Cutoff", 1e-6, "MaxDim", 1));
		psi_defined = true;
		//Compute the density matrix rho associated to the pure state |psi>
		C.psi2rho(psi, argsRho);
		cout << "psi2rho done.\n";
		cout.flush();
	}

	Cplx tr = C.trace_rho();
	cout << "Tr{rho} before re-normalization: " << tr << endl;
	C.rho /= tr; //Normalize rho so that Tr[rho]=1 (otherwise we would have Tr[rho^2]=MPS norm=1)
	tr = C.trace_rho();
	cout << "Tr{rho} after re-normalization: " << tr << endl;

	cout << "Tr{rho^2} =";
	cout.flush();
	const Cplx tr2 = C.trace_rho2();
	cout << tr2 << endl;
	if (psi_defined)
	{
		if (std::abs(tr - 1) > 1e-1 || std::abs(tr2 - 1) > 1e-1)
			cerr << "Error, these traces should be 1 for a pure state |psi><psi|.\n", C.rho /= tr;
		//Check a few simple observables, using rho and psi
		vector<string> ops = {"Sz", "S+", "S-", "Sx", "Sy"};
		vector<double> psi_factor = {2., 1., 1., 2., 2.};

		double err = 0;
		for (int i = 1; i <= N; i++)
		{
			for (unsigned int o = 0; o < ops.size(); o++)
			{
				const string &opname = ops[o];
				Cplx with_rho = C.Expect(opname, i);
				Cplx with_psi = psi_factor[o] * PureStateObs(opname, psi, i, C.sites);
				err += std::abs(with_rho - with_psi);
				if (std::abs(with_rho - with_psi) > 1e-2)
				cerr << "Error: <psi|" << opname << "(" << i << ")|psi>=" << with_psi << "\t"
					<< "Tr[rho*" << opname << "(" << i << ")]=" << with_rho << endl,
					exit(1);
			}
		}
	//err /= N * ops.size();
	//cout << "Compare observables " << ops << " in |psi> and rho: average precision=" << err << endl;
	}
	//-----------------------------------------------------
	//Construct the Lindbladian from the parameters (unitary and dissipative terms)

	SetLindbladian(C, param, lattice);

	//-----------------------------------------------------
	// Compute Lindblad^dagger * Lindblad
	/*
	if (param.longval("dmrgLDL") == 1)
	{
	MPO LD = toMPO(C.LindbladianDag);
	MPO L = toMPO(C.Lindbladian);
	//MPO LDL=toMPO(C.Lindbladian);nmultMPO(L, LD, LDL); // Now LDL=L^Dager * L    itensor v2

	MPO LDL = nmultMPO(L, prime(LD));
	LDL.mapPrime(2, 1); //itensor v3

	cout << "Maximum bond dimension of L^dag*L (MPO):" << maxLinkDim(LDL) << endl;
	auto sweeps = Sweeps(param.longval("sweepsRho"));
	//Specify max number of states kept each sweep
	const int m = param.longval("max_dim_rho");
	if (param.stringval("load_state_file") != "")
	{ //if we are starting from some previous rho, start with the largest allowed bond dimension
	  sweeps.maxdim() = m;
	}
	else
	{
	  //otherwise, increase gradually the matrix dimensions dunring the first sweeps
	  sweeps.maxdim() =
	      min(5, m),
	  min(5, m), min(10, m), min(10, m), min(20, m), min(20, m),
	  min(50, m), min(50, m), min(50, m),
	  min(100, m), min(100, m), min(100, m),
	  min(150, m), min(150, m), min(150, m),
	  min(200, m), min(200, m), min(200, m),
	  min(300, m), min(300, m), min(300, m),
	  min(400, m), min(400, m), min(400, m),
	  min(500, m), min(500, m), min(500, m), m;
	}
	sweeps.cutoff() = param.val("cut_off_rho");
	//Run the DMRG
	//Convergence criterium on the energy passed to the DMRGObserver.
	//3rd argument=true in the constructor of MyDMRGObserver => relative variations |dE/E| are considered
	MyDMRGObserver obs(C.rho, param.val("LDLconv"), true);
	dmrg(C.rho, LDL, sweeps, obs, "Quiet");
	Cplx z = C.trace_rho();
	C.rho /= z; //Normalize rho so that Tr[rho]=1 (otherwise we would have Tr[rho^2]=MPS norm=1)
	}
	*/
	//-----------------------------------------------------
	cout << "Compute exp(i*tau*Linblad) as an MPO... ";
	cout.flush();

	const double tau = param.val("tau");
	const int o = param.val("Trotter_order");
	TimeEvolver TE; //Object defined in "TimeEvolution.h" and "TimeEvolution.cc"
	TE.init(tau, C.Lindbladian, argsRho, o);
	cout << "done.\n";
	cout.flush();
	cout << "Largest bond dimension of exp(tau*L)  (MPO):" << maxLinkDim(TE.expL1) << endl;
	double ttotal = param.val("t_final");
	const int nt = int(ttotal / tau + (1e-9 * (ttotal / tau)));

	string output_prefix = param.stringval("output_files_prefix");
	output_prefix += ".N=" + to_string(N);
 
	//-----------------------------------------------------
	// Some preparation/checks for the 1-qbit observables
	ofstream file_1q(output_prefix + ".obs.1q.dat");
	file_1q << "#Component\ttime_t\tsite_i\tExpectationValue" << endl;
	auto components = param.stringvec("1q_components");
	for (auto &s : components)
	{
		if (s.length() != 1)
			cerr << "Error: " << s << " is an unknown 1-qubit component (should be in {x,y,z} or in {X,Y,Z}).\n", exit(1);
		char c = toupper(s[0]);
		if (c != 'X' && c != 'Y' && c != 'Z')
			cerr << "Error: " << s << " is an unknown 1-qubit component (should be in {x,y,z} or in {X,Y,Z}).\n", exit(1);
	}
	vector<long> sit = param.longvec("1q_sites");
	if (sit.size() == 0)
	{ //If no sites are given explicitely we consider all: 1,...,N
	    sit.resize(N);
		iota(sit.begin(), sit.end(), 1);
	}
	for (int i : sit)
	{
		if (i < 1 || i > N)
			cerr << "Error: invalid site i=" << i << " found in list `1q_sites`.\n", exit(1);
	}
	file_1q.precision(15);
	//-----------------------------------------------------
	// Some preparation/checks for the 2-qbit observables
	ofstream file_2q(output_prefix + ".obs.2q.dat");
	file_2q << "#Component\ttime_t\tsite_i\tsite_j\tExpectationValue" << endl;
	file_2q.precision(15);
	vector<long> sit2 = param.longvec("2q_sites");
	if (sit2.size() % 2 == 1)
		cerr << "Error: the list of sites given in the parameter `2q_sites` should have an even length.\n", exit(1);

	if (sit2.size() == 0)
	{ //If no sites are given explicitely we consider all pairs 1,2,1,3,...,1,N,    2,1,2,3,2,4,...,2,N,  ...  N,N-1
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
			{
				if (i != j)
					sit2.push_back(i), sit2.push_back(j);
			}
	}
	for (unsigned int n = 0; n < sit2.size(); n += 2)
	{
		const int i = sit2[n], j = sit2[n + 1];
		if (i < 1 || i > N)
			cerr << "Error: invalid site i=" << i << " found in list `2q_sites`.\n", exit(1);
		if (j < 1 || j > N)
			cerr << "Error: invalid site i=" << i << " found in list `2q_sites`.\n", exit(1);
	}

	auto components2 = param.stringvec("2q_components");
	for (auto &s : components2)
	{
		if (s.length() != 2)
			cerr << "Error: " << s << " is an unknown 2-qubit component (should be a pair in (x,y,z)*(x,y,z)).\n", exit(1);
		for (int n = 0; n <= 1; n++)
		{
			char c = toupper(s[n]);
			if (c != 'X' && c != 'Y' && c != 'Z')
				cerr << "Error: " << s << " is an unknown component (should be a pair in (x,y,z)*(x,y,z)).\n", exit(1);
		}
	}
	//-----------------------------------------------------
	ofstream file_ent(output_prefix + ".entropy.dat");
	file_ent << "#time\tS_2\tOSEE(center)\tdBondDim(center)\tdBondDim(max)" << endl;
	file_ent.precision(15);
	//-----------------------------------------------------
	const int obs = param.longval("output_step");
	auto start_evolve = steady_clock::now();
	auto duration = duration_cast<milliseconds>(start_evolve - start_sim);
	cout << "\nSimulation duration so far: " << duration.count() / 1000. << "s" << endl;
	auto prev_step = steady_clock::now();

	char buf[100];
	for (int n = 0; n <= nt; n++)
	{
		if (obs > 0)
		{
			if ((n % obs) == 0 || n == nt)
			{
				auto time_step = steady_clock::now();
				duration = duration_cast<milliseconds>(time_step - prev_step);
				cout << "\nTime since previous printout: " << duration.count() / 1000. << "s";
				auto tot_duration = duration_cast<seconds>(time_step - start_sim);
				sprintf(buf, "%.2fhr", tot_duration.count() / 3600.);
				cout << ".\tTotal simulation duration: " << buf << endl;
				prev_step = time_step;

				const double t = n * tau;
				//Some data about this time step
				cout << "Solution time t=" << n * tau << "\t----------------------------\n";

				if (param.boolval("b_force_rho_Hermitian") != 0)
					C.MakeRhoHermitian(argsRho);

				const Cplx tr2 = C.trace_rho2();
				const Real osee = OSEE(C.rho, N / 2); //This is not the entanglement entropy, but the operator-space entanglement entropy (OSEE) of rho, associated to a cut in the middle of te system (N/2)

				const double S_2 = 1.0 / (1.0 - 2.0) * log(tr2.real());
				const int bd = BondDim(C.rho, N / 2), bd_max = maxLinkDim(C.rho);
//				cout << "\tTr{rho}=" << C.trace_rho() << "\tTr{rho^2} =" << tr2 << "\tRényi Entropy S_2 =" << S_2
				cout << "\tTr{rho}=" << C.trace_rho() << "\tRényi Entropy S_2 =" << S_2
				     << "\n\tCenter bond dimension: " << bd << "\tMax bond dimension of rho: " << bd_max
				     << "\n\tOperator Space Entropy at center bond :" << osee << endl;
				file_ent << t << " \t" << S_2 << " \t" << osee << " \t" << bd << " \t" << bd_max << endl;
				//-----------------------------------------------------------------------------------------
				//compute the 1-q  observables and write them to file_1q
				cout << "\tObservables:";
				int count = 0;

				for (long &i : sit)
				{
					for (auto &s : components)
					{
					  Cplx expectation_value;
					  if (tolower(s[0]) == 'x')
					    expectation_value = C.Expect("Sx", i);
					  if (tolower(s[0]) == 'y')
					    expectation_value = C.Expect("Sy", i);
					  if (tolower(s[0]) == 'z')
					    expectation_value = C.Expect("Sz", i);
					  if (abs(expectation_value.imag()) > 1e-3)
					    cout << "Warning: <S^" << s << "(" << i << ")>=" << expectation_value << " is not real\n";
					  file_1q << char(toupper(s[0])) << "\t" << t << "\t" << i
					          << "\t" << expectation_value.real() << endl;
					  count++;
					}
//					file_1q << endl;
				}
				cout << "\n\t\t" << count << " 1-qubit expectation values have been computed and written to a file.";
				file_1q << endl; //Skip a line between each time step
				file_1q.flush();
		        //-----------------------------------------------------------------------------------------
				//compute the 2-q  observables and write them to file_2q
				count = 0;
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
						if (abs(expectation_value.imag()) > 1e-3)
						cout << "Warning: <" << c1 << "(" << i << ")" << c2 << "(" << j << ")>=" << expectation_value << " is not real.\n";
						file_2q << char(toupper(s[0])) << char(toupper(s[1])) << "\t" << t << "\t" << i << "\t" << j << "\t" << expectation_value.real() << endl;
						count++;
					}
				}
				cout << "\n\t\t" << count << " 2-qubit expectation values (correlations) have been computed and written to a file.\n";
				file_2q << endl; //Skip a line between each time step
				file_2q.flush();
	        //-----------------------------------------------------------------------------------------
			}
			else
			{
				cout << ".";
				cout.flush();
			}
		}
		if (n < nt)
		{
			TE.evolve(C.rho);
			Cplx z = C.trace_rho(); //Should be very close to 1, since the Lindblad evolution preserves Tr[rho]
			if (std::abs(z - 1) > 1e-2)
				cout << "Warning: Tr[rho]<>1 :" << z << endl;
			if (param.val("b_force_rho_trace") != 0)
				C.rho /= z;
		}
	}

	bool b_save_state = param.boolval("b_save_final_state");
	if (b_save_state)
	{
		string f1 = output_prefix + ".state.ops";
		writeToFile(f1, C.siteops);
		string f2 = output_prefix + ".state.rho";
		writeToFile(f2, C.rho);
		string f3 = output_prefix + ".state.sites";
		writeToFile(f3, C.sites);
		cout << "The final state was saved to disk, using 3 files:\n" << f1 << endl << f2 << endl << f3 << endl;
	}
	cout << endl;
	auto end_sim = steady_clock::now();
	auto tot_duration = duration_cast<seconds>(end_sim - start_sim);
	sprintf(buf, "%.2fhr", tot_duration.count() / 3600.);
	cout << "\nTotal simulation duration: " << buf << endl;
	return 0;
}

// Old initialization code
/*
  if (param.stringval("load_purestate_file") != "" && param.stringval("load_state_file") != "")
    cerr << "Error, conflict in parameters:load_purestate_file=" << param.stringval("load_purestate_file")
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
          { // Start the DMRG from the Néel state (|Néel> = |up|down|up|down|up|...> )
            if (i % 2 == 1)
              initState.set(i, "Dn");
            else
              initState.set(i, "Up");
          }
        }
        else
        {
          if (param.longval("up_init") != 0 && param.longval("down_init") != 0)
            cerr << "Error: conflicting initialization options for the DMRG:up_init and  down_init.\n", exit(1);

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
        cout << "Initial energy=" << energy << endl;
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
          cout << "the final pure state was written to disk, in files " << f1 << ", " << f2 << " and " << f3 << ".\n";
        }
      }
      else
      { //read psi from a file
        string file_name = param.stringval("load_purestate_file");
        file_name += "_N=" + to_string(N) + ".psi";
        psi = MPS(C.sites);
        cout << "Read the initial pure state (wave-function) from file '" << file_name << "'...";
        cout.flush();
        readFromFile(file_name, psi);
        cout << "done.\n";
      }
      psi_defined = true;

      //Compute the density matrix rho associated to the pure state |psi>
      C.psi2rho(psi, argsRho);
      cout << "psi2rho done.\n";
      cout.flush();
    }
    else
    { //Read the density matrix from disk

      string file_name = param.stringval("load_state_file");
      file_name += "_N=" + to_string(N) + ".rho";
      cout << "Read the initial rho from the file '" << file_name << "'...";
      cout.flush();
      readFromFile(file_name, C.rho);
      cout << "done.\n";
      if (param.val("InitialOrthoRho") != 0)
      {
        cout << "C.rho.orthogonalize...";
        cout.flush();
        C.rho.orthogonalize(Args("Cutoff", param.val("cut_off_rho"), "MaxDim", param.longval("max_dim_rho")));
        cout << "done.\n";
        cout.flush();
      }
    }
  }
  else
  {
    //Start from rho ~ identity (infinite temperature density matrix)
    cout << "Initialize rho ~ identity (infinite temperature density matrix).\n";
    C.rho = C.Identity;
  }
*/