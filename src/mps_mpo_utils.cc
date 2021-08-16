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

#include "mps_mpo_utils.h"

//____________________________________________________________________
double Entropy(MPS psi, int i) //returns the von Neumann entropy on some bond (i,i+1)
{
  auto bond_index=commonIndex(psi.A(i),psi.A(i+1),"Link");
  int bond_dim=bond_index.dim();
  
  psi.position(i); 
  ITensor wf = psi.A(i)*psi.A(i+1);
  auto U = psi.A(i);
  ITensor S,V;
  //Remark: We know that the rank of wf is at most bond_dim, so we specify
  //this value to the SVD routine, in order to avoid many spurious small singular values (like ~ 1e-30)
  //which should in fact be exaclty zero.
  auto spectrum = svd(wf,U,S,V,{"MaxDim",bond_dim});
  Real SvN = 0.;
  Real sum=0;
//  cout<<"\tSingular value decomposition: dim="<<spectrum.numEigsKept()<<endl;
//  cout<<"\t\tLargest sing. val:"<<spectrum.eig(1);
//  cout<<",\tsmallest sing. val:"<<spectrum.eig(spectrum.numEigsKept());
  for(auto p : spectrum.eigs()) {
    sum+=p;
    SvN += -p*log(p);
  }
//  cout<<",\tsum ="<<sum<<endl;
  return SvN;
}
//____________________________________________________________________
double OSEE(MPS rho, int i) {//rho is a density matrix (in MPS form)
  const Cplx tr2=innerC(rho,rho);
  const double s= Entropy(rho,i);
  return s/tr2.real()+log(tr2.real());
}
//____________________________________________________________________
void prints_SVD_spectrum(ostream& o,MPS psi, int i) {
  auto bond_index=commonIndex(psi.A(i),psi.A(i+1),"Link");
  int bond_dim=bond_index.dim();
  
  psi.position(i); 
  ITensor wf = psi.A(i)*psi.A(i+1);
  auto U = psi.A(i);
  ITensor S,V;
  auto spectrum = svd(wf,U,S,V,{"MaxDim",bond_dim});
  double sum=0;
  for(auto p : spectrum.eigs()) sum+=p;
  int dim=spectrum.numEigsKept();
  o<<"dim="<<dim<<"\tsum ="<<sum<<endl;
  double sum2=0,p0=0;int n=0;bool middle=false;
  for(auto p : spectrum.eigs()) {
    sum2+=p;
    if (n==0) p0=p;
    if (n<5 || (dim-n)<5) {
      o<<"p["<<n<<"]="<<p<<"\tmissing weight:"<<double(1.0)-sum2/sum<<"\tp["<<n<<"]/p[0]="<<p/p0<<endl;
    } else {
      if (!middle) o<<"..."<<endl,middle=true;
    }
    n++;
  }
  
}
//____________________________________________________________________
