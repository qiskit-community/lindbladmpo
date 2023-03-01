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

#ifndef _PAULI_
#define _PAULI_

#include "itensor/all.h"
#include <string>

using namespace itensor;
using namespace std;


//____________________________________________________________________

class PauliSite {
  // Index which can take 4 values
  // These values correspond to some natural basis of a  4-dimensionnal local Hilbert space
  // spanned by the 3 Pauli matrices, and the Identity matrix.
  // This space is used to encode the density matrix of a spin-1/2 model.
  Index s;
 public:
  PauliSite();
  PauliSite(Index I);
  PauliSite(Args const& args = Args::global());
  Index index() const;
  IndexVal state(string const& state) const;
  ITensor op(const string& opname, const Args&args) const;
};

using Pauli = BasicSiteSet<PauliSite>;

//____________________________________________________________________
class SpinHalfSystem {
  
 public:
  const int N;//number of sites;

  SpinHalf sites;//SiteSet for pure-state descrition
  
  //siteop: chain of 2*2 matrices (on-site operators)
  //a wave-function for this "doubled" chain represents a density matrix

  Pauli siteops; //SiteSet for mixed-state descrition (density matrix rho)
  
  //Identity (density matrix) operator (=wave function for the doubled-chain)
  MPS Identity;
  
  //Density matrix
  MPS rho;

  //Lindbladian operator (coded as an autoMPO)
  AutoMPO Lindbladian;

  //Adds a single-spin term to the Lindbladian above
  void AddSingleSpinBath(double GammaPlus, double GammaMinus, double GammaDephasing, int site);
  
  //Trace of rho
  complex<double> trace_rho() const;
  complex<double> trace_rho2() const;

  //Constructor
  SpinHalfSystem(int );//argument = number of spins

  //Initialize the density matrix
  void ConstructIdentity();

  //Convert a pure state (psi) of the chain into a density matrix rho (projector onto psi)
  //this is usefull to build the density matrix rho at t=0 if we want
  //to start the simulation from the ground-state |psi> of some Hamiltonian.
  void psi2rho(const MPS& psi,const Args& args = Args::global());

  //Expectation value of some single-site operator
  Cplx Expect(const string& opname,int i) const;
  //Expectation value of some two-site operator
  Cplx Expect(const string& opname1,int i1,const string& opname2,int i2) const;
  //Expectation value of some three-site operator
  Cplx Expect(const string& opname1,int i1,const string& opname2,int i2,const string& opname3,int i3) const;
  
  //Expectation value of some product of Pauli operators
  Cplx Expect(const vector <string>& opnames,const vector <int>& indices) const;
  
  void MakeRhoHermitian(Args args = Args::global()); 

};
#endif
