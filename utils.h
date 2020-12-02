#ifndef _UTILS_
#define _UTILS_
#include "itensor/all.h"
#include <string>
#include <vector>
#include <map>
using namespace itensor;
using namespace std;
//____________________________________________________________
//The function below translate numbers (etc.) into character strings
//the second parameter (optional) is the precision (digits)

template <class T>
inline string to_string(const T &t, unsigned int precision = 0)
{
  stringstream ss;
  if (precision > 0)
    ss.precision(precision);
  ss << t;
  return ss.str();
}
//____________________________________________________________
template <class T>
ostream &operator<<(ostream &o, const vector<T> &v)
{
  for (unsigned int i = 0; i < v.size(); i++)
    o << "[" << i << "]" << v[i] << " ";
  return o;
}
//____________________________________________________________

inline double char2double(char *a)
{
  char *end_ptr;
  const double x = strtod(a, &end_ptr);
  if (end_ptr == a || ('\0' != *end_ptr))
    cout << endl
         << "ERROR :" << a
         << " is not a valid format for a double."
         << endl,
        exit(0);
  return x;
}
//____________________________________________________________
class Parameters : public map<string, double>
{
public:
  double val(string var_name) const
  {
    map<string, double>::const_iterator it = find(var_name);
    if (it == end())
    {
      cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
      return 0;
    }
    else
      return it->second;
  }
  long longval(string var_name) const
  {
    double v = val(var_name);
    if (abs(double(round(v)) - v) < 1e-6)
    {
      return long(round(v));
    }
    else
    {
      cout << "Error, parameter " << var_name << "=" << v << " is not a long" << endl, exit(0);
      return 0;
    }
  }
  void PRint(ostream &o) const
  {
    for (map<string, double>::const_iterator it = begin(); it != end(); it++)
    {
      o << it->first << "=" << it->second << endl;
    }
  }
  void ReadArguments(int argc, char *argv[])
  {
    for (int n = 1; n < argc; n++)
    {
      string var_name(argv[n]);
      map<string, double>::const_iterator it = find(var_name);

      if (it != end())
      {
        n++;
        if (n == argc)
          cerr << "Error: missing value after " << var_name << endl, exit(0);
        operator[](var_name) = char2double(argv[n]);
      }
      else
      {
        cerr << "Syntax error :" << var_name << endl;
        cout << "List of command-line parameters :\n";
        PRint(cout);
        exit(0);
      }
    }
  }
};
//____________________________________________________________

class MyDMRGObserver : public DMRGObserver
{
  double previous_energy;
  const double precision;
  bool relative; // false => stop criterium on |E_n - E_{n-1} |, true => criterium on |E_n - E_{n-1} |/|E_n|
public:
  MyDMRGObserver(const MPS &psi, double prec = 1e-10, bool rel = false) : DMRGObserver(psi), precision(prec), relative(rel) {}

  bool checkDone(const Args &args = Args::global())
  {
    const double energy = args.getReal("Energy", 0);
    cout << "    Change in <L^+ L>: " << energy - previous_energy << "\trelative change: " << (energy - previous_energy) / (abs(previous_energy) + 1e-15) << endl;
    if (relative)
    {
      if (abs(energy - previous_energy) / (abs(previous_energy) + 1e-15) < precision)
      {
        cout << "   <L^+ L> has converged (relative change <" << precision << ") => stop the DMRG.\n";
        return true;
      }
      else
      {
        previous_energy = energy;
        return false;
      }
    }
    else
    {
      if (abs(energy - previous_energy) < precision)
      {
        cout << "   Energy has converged (change <" << precision << ") => stop the DMRG.\n";
        return true;
      }
      else
      {
        previous_energy = energy;
        return false;
      }
    }
  }
};
//_____________________________________________________
inline Cplx PureStateObs(const string &opname, MPS psi, int i, const SiteSet &sites) //<psi|O|psi> at some site i
{
  auto O = op(sites, opname, i);
  psi.position(i);
  ITensor A = psi.A(i);
  ITensor A_ = dag(prime(psi.A(i), "Site"));
  A *= O;
  ITensor B = A * A_;
  return B.cplx();
}
//____________________________________________________________________

inline int BondDim(const MPS &psi, int i)
{
  auto bond_index = commonIndex(psi.A(i), psi.A(i + 1), "Link");
  return (bond_index.dim());
}
inline int BondDim(const MPO &H, int i)
{
  auto bond_index = commonIndex(H.A(i), H.A(i + 1), "Link");
  return (bond_index.dim());
}
//____________________________________________________________________
double Entropy(MPS psi, int i); //returns the von Neumann entropy on some bond (i,i+1)
//____________________________________________________________________
double OSEE(MPS rho, int i);
void prints_SVD_spectrum(ostream &o, MPS psi, int i);
//____________________________________________________________________
#endif
