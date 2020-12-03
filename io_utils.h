#ifndef _IO_UTILS_
#define _IO_UTILS_
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

#endif
