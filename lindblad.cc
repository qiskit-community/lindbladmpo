#include "itensor/all.h"

#include "Pauli.h"
#include "io_utils.h"
#include "mps_mpo_utils.h"
#include "SimpleSquareLattice.h"
#include "TimeEvolution.h"
#include "ModelParameters.h"
#include "MyLindbladian.h"
using namespace itensor;
using namespace std;

#define UNDEFINED_VAL 9999

//_____________________________________________________
int main(int argc, char *argv[])
{

  ModelParameters param;

  //Read some input parameters from the command line (can be empty -> default values for all parameters):
  param.ReadArguments(argc, argv);

  //If a filename is given as the parameter "inputfile" => some input parameters are read from that file
  string inputfilename = param.stringval("inputfile");
  if (inputfilename != "")
    param.ReadFromFile(inputfilename);

  //Now param contains all the the parameters, default values (see SimulationParameters.h and ModelParameters.h) or those provided on the command-line
  cout << endl;
  cout.precision(10);
  param.check();
  param.PRint(cout);

  int Lx = param.longval("Lx");
  int Ly = param.longval("Ly");
  const int N = Lx * Ly;
  Lattice2d lattice(Lx, Ly, param.boolval("b_periodic_x"), param.boolval("b_periodic_y"));
  SpinHalfSystem C(N);
  if (param.longval("read_wf") != 0 || param.longval("read_rho") != 0)
  {
    readFromFile("sites_file", C.sites);
    readFromFile("siteops_file", C.siteops);
    cout << "C.sites and C.siteops were read from 'sites_file' and 'siteops_file'.\n";
    C.rho = MPS(C.siteops);
    C.Lindbladian = AutoMPO(C.siteops);
    C.LindbladianDag = AutoMPO(C.siteops);
  }
  C.ConstructIdentity(); //Construct the density matrix corresponding to inifinite temperature (~ Identity)

  auto args = Args("Cutoff", param.val("Cutoff"), "MaxDim", param.longval("MaxDim"));

  //Note on the option below: if we write "Normalize=true" (default for fitapplyMPO)
  //we would normalize rho such that Tr[rho^2]=1, (norm of the MPS)
  //which is of course not appropriate for a density matrix.
  //auto argsRho=args;
  //argsRho.add("Normalize",false);argsRho.add("MaxDim",param.longval("MaxDimRho"));argsRho.add("Cutoff",param.longval("CutoffRho"));

  auto argsRho = args;
  argsRho.add("Normalize", false);
  argsRho.add("MaxDim", param.longval("MaxDimRho"));
  argsRho.add("Cutoff", param.val("CutoffRho"));

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

  //-----------------------------------------------------
  MPS psi;
  bool psi_defined = false;

  if (param.longval("read_wf") != 0 && param.longval("read_rho") != 0)
    cerr << "Error: conflict in parameters (read_wf+read_rho): should I read psi or rho from the disk?\n", exit(1);

  if (param.longval("rho_inf_init") == 0)
  {
    if (param.longval("read_rho") == 0)

    { //Start from a wave-function (pure state)
      if (param.longval("read_wf") == 0)
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
        const int m = param.longval("MaxDim");
        //We specify bvelow how the max. bond dimension should be increased along the DMRG sweeps
        sweeps.maxdim() = min(5, m), min(5, m), min(10, m), min(10, m), min(10, m), min(20, m), min(20, m), min(50, m), min(50, m), min(100, m), min(100, m), min(200, m), min(200, m), m;
        sweeps.cutoff() = param.val("Cutoff");
        //Run the DMRG algorithm
        MyDMRGObserver obs(psi, param.val("energy")); //Convergence criterium on the energy passed to the DMRGObserver

        const double energy = dmrg(psi, H0, sweeps, obs, "Quiet");
        cout << "Initial energy=" << energy << endl;
        if (param.val("write_wf") != 0)
        {
          writeToFile("sites_file", C.sites);
          writeToFile("siteops_file", C.siteops);
          writeToFile("psi_file", psi);
          cout << "the final wave-function was written to disk.\n";
        }
      }
      else
      { //read psi from a file
        cout << "Read the initial wave-function from file...";
        psi = MPS(C.sites);
        readFromFile("psi_file", psi);
      }
      psi_defined = true;

      //Compute the density matrix rho associated to the pure state |psi>
      C.psi2rho(psi, argsRho);
      cout << "psi2rho done.\n";
      cout.flush();
    }
    else
    { //Read the density matrix from disk
      cout << "Read the initial rho from file...";
      cout.flush();
      readFromFile("rho_file", C.rho);
      cout << "done.\n";
      if (param.val("InitialOrthoRho") != 0)
      {
        cout << "C.rho.orthogonalize...";
        cout.flush();
        C.rho.orthogonalize(Args("Cutoff", param.val("CutoffRho"), "MaxDim", param.longval("MaxDimRho")));
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
    if (fabs(tr - 1) > 1e-1 || fabs(tr2 - 1) > 1e-1)
      cerr << "Error, these traces should be 1 for a pure state |psi><psi|.\n", C.rho /= tr;
    //Check a few simple observables, using rho and psi
    vector<string> ops = {"Sz", "S+", "S-", "Sx", "Sy"};
    double err = 0;
    for (int i = 1; i <= N; i++)
    {
      for (unsigned int o = 0; o < ops.size(); o++)
      {
        const string &opname = ops[o];
        Cplx with_rho = C.Expect(opname, i);
        Cplx with_psi = PureStateObs(opname, psi, i, C.sites);
        err += fabs(with_rho - with_psi);
        if (fabs(with_rho - with_psi) > 1e-2)
          cerr << "Error: <psi|" << opname << "(" << i << ")|psi>=" << with_psi << "\t"
               << "Tr[rho*" << opname << "(" << i << ")]=" << with_rho << endl,
              exit(1);
      }
    }
    err /= N * ops.size();
    cout << "Compare observables " << ops << " in |psi> and rho: average precision=" << err << endl;
  }
  //-----------------------------------------------------
  //Construct the Lindbladian from the parameters (unitary and dissipative terms)

  SetLindbladian(C, param, lattice);

  //-----------------------------------------------------
  // Compute Lindblad^dagger * Lindblad
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
    const int m = param.longval("MaxDimRho");
    if (param.longval("read_rho") != 0)
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
    sweeps.cutoff() = param.val("CutoffRho");
    //Run the DMRG
    //Convergence criterium on the energy passed to the DMRGObserver.
    //3rd argument=true in the constructor of MyDMRGObserver => relative variations |dE/E| are considered
    MyDMRGObserver obs(C.rho, param.val("LDLconv"), true);
    dmrg(C.rho, LDL, sweeps, obs, "Quiet");
    Cplx z = C.trace_rho();
    C.rho /= z; //Normalize rho so that Tr[rho]=1 (otherwise we would have Tr[rho^2]=MPS norm=1)
  }
  //-----------------------------------------------------
  cout << "Compute exp(i*tau*Linblad) as an MPO... ";
  cout.flush();

  const double tau = param.val("tau");
  const int o = param.val("TrotterOrder");
  TimeEvolver TE; //Object defined in "TimeEvolution.h" and "TimeEvolution.cc"
  TE.init(tau, C.Lindbladian, argsRho, o);
  cout << "done.\n";
  cout.flush();
  cout << "Largest bond dimension of exp(tau*L)  (MPO):" << maxLinkDim(TE.expL1) << endl;
  double ttotal = param.val("T");
  const int nt = int(ttotal / tau + (1e-9 * (ttotal / tau)));

  string name = param.FileName();

  ofstream file_sz("sz_profile." + name); //This file will the full { <S^z(i)> i=1...N} profile as a function of time
  file_sz << "#\"t=time\"" << endl;
  file_sz << "#time\t1\t<S^z(1)>\t<S^x(1)>" << endl;
  file_sz << "#time\ti\t<S^z(i)>\t<S^x(i)>" << endl;
  file_sz << "#time\tN\t<S^z(N)>\t<S^x(N)>" << endl
          << endl;
  file_sz.precision(15);

  ofstream file_m("time_series." + name); //This file will contain various observables (in different columns) as a function of time
  file_m << "#time\t<S^z>\t<S^x>\t<S^y>\tS_2\tOSEE(center)\tdBondDim(center)\tdBondDim(max)"
         << "\t\\eta_{zz}\t\\eta_{xx}\t\\eta_{yy}\t\\eta_{xz}\t\\eta_{yz}\t\\eta_{xy}" << endl;
  file_m.precision(15);

  const int obs = param.longval("obs");
  for (int n = 0; n <= nt; n++)
  {
    if (obs > 0)
    {
      if ((n % obs) == 0 || n == nt)
      {
        const double t = n * tau;
        //Some data about this time step
        cout << "\nt=" << n * tau << "\t----------------------------\n";

        if (param.longval("b_force_rho_Hermitian") != 0)
          C.MakeRhoHermitian(argsRho);

        const Cplx tr2 = C.trace_rho2();
        const Real osee = OSEE(C.rho, N / 2); //This is not the enanglement entropy, but the operator-space entanglement entropy (OSEE) of rho, associated to a cut in the middle of te system (N/2)

        const double S_2 = 1.0 / (1.0 - 2.0) * log(tr2.real());
        const int bd = BondDim(C.rho, N / 2), bd_max = maxLinkDim(C.rho);
        cout << "\tTrace[rho]=" << C.trace_rho() << "\tTr{rho^2} =" << tr2 << "\tRényi Entropy S_2 =" << S_2
             << "\tCenter-bond dimension for rho:" << bd << "\tMax bond dimension of rho:" << bd_max
             << "\n\tOperator Space Entropy @center bond :" << osee;

        cout << "\n\tMagnetization:\n\t\t";
        {
          file_sz << "\"t=" << t << "\"\n";
          for (int i = 1; i <= N; i++)
          {
            const Cplx sz = C.Expect("Sz", i);
            if (abs(sz.imag()) > 1e-3)
              cout << "Warning: <Sz(" << i << ")>=" << sz << " is not real\n";
            const Cplx sx = C.Expect("Sx", i);
            if ((i == 1) || abs(i - N / 2) < 1 || i == N)
              cout << "Sz(" << i << ")=" << sz << "\t";
            file_sz << t << "\t" << i - N / 2 << "\t" << sz.real() << "\t" << sx.real() << endl;
          }
        }
        double sz_tot = 0, sz_center = 0, sza = 0, szb = 0, sxa = 0, sxb = 0, sya = 0, syb = 0,
               szsz = 0, sxsx = 0, sysy = 0, sysz = 0, sxsz = 0, sxsy = 0;
        const int center_site = N / 2;
        int a = N / 2, b = a + 1;
        if (a < 1)
          a = 1;
        if (b > N)
          b = N;
        {
          for (int i = 1; i <= N; i++)
          {
            const double s = C.Expect("Sz", i).real();
            sz_tot += s;
            if (i == center_site)
              sz_center += s;
            if (i == a)
              sza = s, sxa = C.Expect("Sx", i).real(), sya = C.Expect("Sy", i).real();
            if (i == b)
              szb = s, sxb = C.Expect("Sx", i).real(), syb = C.Expect("Sy", i).real();
          }
          cout << "<Sz_center>=" << sz_center << "\t"
               << "<Sz_tot>=" << sz_tot
               << "\n\t\t<n_center>=" << sz_center + 0.5 << "\t"
               << "<n>=" << sz_tot / N + 0.5
               << "\n\t\t<Sz(a=" << a << ")>=" << sza << "\t<Sz(b=" << b << ")>=" << szb
               << "\n\t\t<Sx(a=" << a << ")>=" << sxa << "\t<Sx(b=" << b << ")>=" << sxb
               << "\n\t\t<Sy(a=" << a << ")>=" << sya << "\t<Sy(b=" << b << ")>=" << syb;
          cout << "\n\tCorrelations:";

          //Two-point correlations
          szsz = C.Expect("Sz", a, "Sz", b).real();
          cout << "\n\t\t<Sz(a)Sz(b)>=" << szsz << "\t<Sz(a)Sz(b)>^c=" << szsz - sza * szb;
          sxsx = C.Expect("Sx", a, "Sx", b).real();
          cout << "\n\t\t<Sx(a)Sx(b)>=" << sxsx << "\t<Sx(a)Sx(b)>^c=" << sxsx - sxa * sxb;
          sysy = C.Expect("Sy", a, "Sy", b).real();
          cout << "\n\t\t<Sy(a)Sy(b)>=" << sysy << "\t<Sy(a)Sy(b)>^c=" << sysy - sya * syb;

          sxsz = C.Expect("Sx", a, "Sz", b).real();
          cout << "\n\t\t<Sx(a)Sz(b)>=" << sxsz << "\t<Sx(a)Sz(b)>^c=" << sxsz - sxa * szb;
          sysz = C.Expect("Sy", a, "Sz", b).real();
          cout << "\n\t\t<Sy(a)Sz(b)>=" << sysz << "\t<Sy(a)Sz(b)>^c=" << sysz - sya * szb;
          sxsy = C.Expect("Sx", a, "Sy", b).real();
          cout << "\n\t\t<Sx(a)Sy(b)>=" << sxsy << "\t<Sx(a)Sy(b)>^c=" << sxsy - sxa * syb << endl;
        }
        file_sz << endl
                << endl;
        file_sz.flush();

        file_m << t << "\t"
               << C.Expect("Sz", center_site).real() << "\t"  //col 2
               << C.Expect("Sx", center_site).real() << "\t"  // 3
               << -C.Expect("Sy", center_site).real() << "\t" //4
               << S_2 << "\t"                                 //5
               << osee << "\t"                                //6
               << bd << "\t"                                  //7
               << bd_max << "\t"                              //8
               << szsz - sza * szb << "\t"                    //9
               << sxsx - sxa * sxb << "\t"                    //10
               << sysy - sya * syb << "\t"                    //11
               << sxsz - sxa * szb << "\t"                    //12
               << -sysz + sya * szb << "\t"                   //13
               << -sxsy + sxa * syb << "\t"                   //14
               << endl;
        file_m.flush();
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
      if (fabs(z - 1) > 1e-2)
        cout << "Warning: Tr[rho]<>1 :" << z << endl;
      if (param.val("b_force_rho_trace") != 0)
        C.rho /= z;
    }
  }

  if (param.val("write_rho") != 0)
  {
    writeToFile("siteops_file", C.siteops);
    writeToFile("rho_file", C.rho);
    writeToFile("sites_file", C.sites);
    cout << "the density matrix rho was written to disk, in files 'rho_file', 'siteops_file' and 'sites_file'.\n";
  }
  cout << endl;
  return 0;
}
