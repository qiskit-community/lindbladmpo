#ifndef _IO_UTILS_
#define _IO_UTILS_
#include "itensor/all.h"
#include <string>
#include <vector>
#include <map>
#include <regex>
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
class Parameters_old : public map<string, double>
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
class Parameters : public map<string, string>
{
public:
  //------------------------------------------------------
  // The method 'ReadFromFile reads an input file (specified by its name) which should be of the following (text) format:
  // variable1 = value_of_variable_1
  // variable2 = value_of_variable_2
  // ...
  // At this stage all values are strings (without any space)
  // and they are added to the Parameters object.
  // Any line without any '=' sign is ignored
  // All the variables should have been previously declared (see SimulationParameters.h or ModelParameters.h)
  // Note: the methods Parameters::val, Parameters::longval, Parameters::boolval, and Parameters::stringval can later be used to retrieve these values as
  // double, long, bool and string.
  // TODO: add support for vectors: vector_variable = x1,x2,x3,x4,x5,...
  void ReadFromFile(string filename)
  {
    ifstream file(filename);
    if (!file) cerr<<"Error: unable to open the file "<<filename<<endl,exit(0);
    else cout<<"Reading parameters from the file "<<filename<<endl;
    string line;
    while (getline(file, line))
    {
      line = regex_replace(line, regex("\\s+"), "");// remove white spaces using a regular expression
      if (line != "")
      {
        string delimiter = "=";
        size_t position = 0;
        string var_name;
        position = line.find(delimiter);
        if (position != string::npos)
        {
          var_name = line.substr(0, position);
          string var_value = line.erase(0, position + delimiter.length());

          map<string, string>::const_iterator it = find(var_name);
          if (it != end())
            operator[](var_name) = var_value;
          else
            cerr << "Error, the input parameter " << var_name << " does not exist in this model.\n", exit(0);
        }
        else
        {
          cout << "Warning: the input line :"
               << line << " has been ignored (no '=')." << endl;
        }
      }
    }
  }
  //------------------------------------------------------
  double val(string var_name) const
  {
    map<string, string>::const_iterator it = find(var_name);
    if (it == end())
    {
      cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
      return 0;
    }
    else
      return stod(it->second);
  }
  //------------------------------------------------------
  long longval(string var_name) const
  {
    map<string, string>::const_iterator it = find(var_name);
    if (it == end())
    {
      cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
      return 0;
    }
    else
    {
      return stoi(it->second);
    }
  }
  //------------------------------------------------------
  long boolval(string var_name) const
  {
    map<string, string>::const_iterator it = find(var_name);
    if (it == end())
    {
      cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
      return 0;
    }
    else
    {
      string s = it->second;
      if (s == "true")
        return true;
      if (s == "TRUE")
        return true;
      if (s == "1")
        return true;
      if (s == "false")
        return false;
      if (s == "FALSE")
        return false;
      if (s == "0")
        return false;
      cout << "Error " << var_name << "=" << it->second << " but a boolean was expected: true/false, TRUE/FALSE or or 1/0\n", exit(0);
    }
  } //------------------------------------------------------
  string stringval(string var_name) const
  {
    map<string, string>::const_iterator it = find(var_name);
    if (it == end())
    {
      cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
      return 0;
    }
    else
      return (it->second);
  }
  //------------------------------------------------------
  void PRint(ostream &o) const
  {
    for (map<string, string>::const_iterator it = begin(); it != end(); it++)
    {
      o << it->first << "=" << it->second << endl;
    }
  }
  //------------------------------------------------------
  void ReadArguments(int argc, char *argv[])
  {
    for (int n = 1; n < argc; n++)
    {
      string var_name(argv[n]);
      map<string, string>::const_iterator it = find(var_name);

      if (it != end())
      {
        n++;
        if (n == argc)
          cerr << "Error: missing value after " << var_name << endl, exit(0);
        operator[](var_name) = string(argv[n]);
      }
      else
      {
        cerr << "Syntax error, the input parameter " << var_name << " does not exist in this model.\n";
        cout << "List of command-line parameters :\n";
        PRint(cout);
        exit(0);
      }
    }
  }
};
#endif
