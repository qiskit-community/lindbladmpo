#include "itensor/all.h"
#include "gates.h"

using namespace itensor;
using namespace std;


//Apply the X gate (qubit i) on a mixed state rho
void ApplyXGate(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"Sx",i);
	rho.ref(i)*=op(siteops,"_Sx",i);
	rho.ref(i).noPrime("Site");
}
//Apply the Y gate (qubit i) on a mixed state rho
void ApplyYGate(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"Sy",i);
	rho.ref(i)*=op(siteops,"_Sy",i);
	rho.ref(i).noPrime("Site");
}
//Apply the X gate (qubit i) on a mixed state rho
void ApplyZGate(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"Sz",i);
	rho.ref(i)*=op(siteops,"_Sz",i);
	rho.ref(i).noPrime("Site");
}
//Apply the SqrtX gate (qubit i) on a mixed state rho
void ApplySqrtXGate(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"SqrtX",i);
	rho.ref(i)*=op(siteops,"_SqrtX",i);
	rho.ref(i).noPrime("Site");
}
void ApplyHGate(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"H",i);
	rho.ref(i)*=op(siteops,"_H",i);
	rho.ref(i).noPrime("Site");
}
void ApplyProjUp(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"projUp",i);
	rho.ref(i)*=op(siteops,"_projUp",i);
	rho.ref(i).noPrime("Site");
}
void ApplyProjDn(MPS &rho, const Pauli &siteops,int i) {
	rho.ref(i)*=op(siteops,"projDn",i);
	rho.ref(i)*=op(siteops,"_projDn",i);
	rho.ref(i).noPrime("Site");
}

void ApplyControlledXYZGateONPureState(MPS &psi, const SpinHalf&sites,int control,int target,string opname,Args args) {
	if (control==target) cerr << "Error, ApplyControlledXYZGateONPureState was called with control=target="<<control<<".\n", exit(1);
	const int i=min(target,control),j=max(target,control);
	const int N=length(psi);	
	{
		MPO gate_mpo(sites);
		//Construct the link indices
		vector<Index> links(N);
		for(int n = 1; n < N; ++n) {
			if (n<i) {//Bond dim. = 1
				links.at(n) = Index(1,format("Link,l=%d",n));
			}
			if (n>=i && n<j) {//Bond dim. = 2
				links.at(n) = Index(2,format("Link,l=%d",n));
			}
			if (n>=j) {//Bond dim. = 1
				links.at(n) = Index(1,format("Link,l=%d",n));
			}
		}
		//Construct 'manually' the tensors of the MPO for the CNOT gate
		for(int n = 1; n <= N; ++n) {
			auto& W = gate_mpo.ref(n);
			if (n==1) {
				Index right=links.at(n);
				W = ITensor(sites(n),prime(sites(n)),right);
				if (n==control) {
					// the bond index of the mpo encodes the information about the satet of the control qubit
					W+=sites.op("projUp",n)  * setElt(right(1)); // |control> = |0> = |up>
					W+=sites.op("projDn",n)  * setElt(right(2)); // |control> = |1> = |dn>
				} else if (n==target) {
					W+=sites.op("Id",n)  * setElt(right(1));// control> = |0>
					W+=2*sites.op(opname,n)  * setElt(right(2));  // Apply 2*"opname" to the target bit if |control> = |1>. The 2* comes from the fact that "Sx" is 0.5*Pauli^x when applied to a Spin-1/2 state (whereas it is equal to 1.0*Pauli^x when applied on au mixed/Pauli/density matrix state)
				} else
				W += sites.op("Id",n) * setElt(right(1)) ;
			}
			if (n>1 && n<N) {
				Index right=links.at(n);
				Index left=links.at(n-1);
				W = ITensor(sites(n),prime(sites(n)),right,left);
				if (n==i && n==control) {
					W+=sites.op("projUp",n)  * setElt(left(1)) * setElt(right(1));
					W+=sites.op("projDn",n)  * setElt(left(1)) * setElt(right(2));
				} else 	if (n==i && n==target) {
					W+=sites.op("Id",n)  * setElt(left(1)) * setElt(right(1));
					W+=2*sites.op(opname,n)  * setElt(left(1)) * setElt(right(2));// Apply 2*opname to the target bit if |control> = |1>
				} else if (n>i && n<j) {
					W+=sites.op("Id",n) * setElt(right(1)) * setElt(left(1)) ;
					W+=sites.op("Id",n) * setElt(right(2)) * setElt(left(2)) ;
				} else if (n==j && n==control) {
					W+=sites.op("projUp",n)  * setElt(left(1)) * setElt(right(1));
					W+=sites.op("projDn",n)  * setElt(left(2)) * setElt(right(1));
				} else if (n==j && n==target) {
					W+=sites.op("Id",n)  * setElt(left(1)) * setElt(right(1));
					W+=2*sites.op(opname,n)  * setElt(left(2)) * setElt(right(1));// Apply 2*opname to target bit if |control> = |1>
				} else
					W+=sites.op("Id",n) * setElt(right(1)) * setElt(left(1)) ;
			}
			if (n==N) {
				Index left=links.at(n-1);
				W = ITensor(sites(n),prime(sites(n)),left);
				if (n==control) {
					W+=sites.op("projUp",n)  * setElt(left(1)) ;
					W+=sites.op("projDn",n)  * setElt(left(2)) ;
				} else if (n==target) {
					W+=sites.op("Id",n)  * setElt(left(1)) ;
					W+=sites.op(opname,n)  * setElt(left(2)) ;// Flip the target bit if |control> = |1>
				} else
					W += sites.op("Id",n) * setElt(left(1)) ;	
			}
		}
		psi=applyMPO(gate_mpo,psi,args);psi.noPrime("Site");
	}
}


void ApplyControlledXYZGate(MPS &rho, const Pauli &siteops,int control,int target,string opname,Args args) {
	if (control==target) cerr << "Error, ApplyControlledXYZGate was called with control=target="<<control<<".\n", exit(1);
	const int i=min(target,control),j=max(target,control);
	const int N=length(rho);	
	for (int braket=0;braket<=1;braket++) {
		MPO gate_mpo(siteops);
		string projUp=(braket==0)?("projUp"):("_projUp");//acting on |ket> or on <bra|
		string projDn=(braket==0)?("projDn"):("_projDn");
		string s=(braket==0)?(opname):("_"+opname);
		//Construct the link indices
		vector<Index> links(N);
		for(int n = 1; n < N; ++n) {
			if (n<i) {//Bond dim. = 1
				links.at(n) = Index(1,format("Link,l=%d",n));
			}
			if (n>=i && n<j) {//Bond dim. = 2
				links.at(n) = Index(2,format("Link,l=%d",n));
			}
			if (n>=j) {//Bond dim. = 1
				links.at(n) = Index(1,format("Link,l=%d",n));
			}
		}
		//Construct 'manually' the tensors of the MPO for the CNOT gate
		for(int n = 1; n <= N; ++n) {
			auto& W = gate_mpo.ref(n);
			if (n==1) {
				Index right=links.at(n);
				W = ITensor(siteops(n),prime(siteops(n)),right);
				if (n==control) {
					// the bond index of the mpo encodes the information about the satet of the control qubit
					W+=siteops.op(projUp,n)  * setElt(right(1)); // |control> = |0> = |up>
					W+=siteops.op(projDn,n)  * setElt(right(2)); // |control> = |1> = |dn>
				} else if (n==target) {
					W+=siteops.op("Id",n)  * setElt(right(1));// control> = |0>
					W+=siteops.op(s,n)  * setElt(right(2));  // Flip the target bit if |control> = |1>
				} else
				W += siteops.op("Id",n) * setElt(right(1)) ;
			}
			if (n>1 && n<N) {
				Index right=links.at(n);
				Index left=links.at(n-1);
				W = ITensor(siteops(n),prime(siteops(n)),right,left);
				if (n==i && n==control) {
					W+=siteops.op(projUp,n)  * setElt(left(1)) * setElt(right(1));
					W+=siteops.op(projDn,n)  * setElt(left(1)) * setElt(right(2));
				} else 	if (n==i && n==target) {
					W+=siteops.op("Id",n)  * setElt(left(1)) * setElt(right(1));
					W+=siteops.op(s,n)  * setElt(left(1)) * setElt(right(2));// Flip the target bit if |control> = |1>
				} else if (n>i && n<j) {
					W+=siteops.op("Id",n) * setElt(right(1)) * setElt(left(1)) ;
					W+=siteops.op("Id",n) * setElt(right(2)) * setElt(left(2)) ;
				} else if (n==j && n==control) {
					W+=siteops.op(projUp,n)  * setElt(left(1)) * setElt(right(1));
					W+=siteops.op(projDn,n)  * setElt(left(2)) * setElt(right(1));
				} else if (n==j && n==target) {
					W+=siteops.op("Id",n)  * setElt(left(1)) * setElt(right(1));
					W+=siteops.op(s,n)  * setElt(left(2)) * setElt(right(1));// Flip the target bit if |control> = |1>
				} else
					W+=siteops.op("Id",n) * setElt(right(1)) * setElt(left(1)) ;
			}
			if (n==N) {
				Index left=links.at(n-1);
				W = ITensor(siteops(n),prime(siteops(n)),left);
				if (n==control) {
					W+=siteops.op(projUp,n)  * setElt(left(1)) ;
					W+=siteops.op(projDn,n)  * setElt(left(2)) ;
				} else if (n==target) {
					W+=siteops.op("Id",n)  * setElt(left(1)) ;
					W+=siteops.op(s,n)  * setElt(left(2)) ;// Flip the target bit if |control> = |1>
				} else
					W += siteops.op("Id",n) * setElt(left(1)) ;	
			}
		}
		rho=applyMPO(gate_mpo,rho,args);rho.noPrime("Site");
	}
}
//Apply the CNOT gate on some mixed state rho, at sites (control,target)
void ApplyCNOTGate(MPS &rho, const Pauli &siteops,int control,int target,Args args) {
	ApplyControlledXYZGate(rho, siteops,control,target,"Sx",args);
}

//Apply the controlled-Z gate on some mixed state rho, at sites (i,j)
void ApplyControlledZGate(MPS &rho, const Pauli &siteops,int i,int j,Args args) {
	ApplyControlledXYZGate(rho, siteops,i,j,"Sz",args);
}

ITensor Hadamard(const SpinHalf& sites, int n) {
    double sqrt05 = pow(.5, .5);
    auto ind = sites(n);
    auto indP = prime(sites(n));
    auto H = ITensor(ind,indP);
    H.set(ind(1),indP(1), sqrt05);
    H.set(ind(1),indP(2), sqrt05);
    H.set(ind(2),indP(1), sqrt05);
    H.set(ind(2),indP(2), -sqrt05);
    return H;       
}

//_____________________________________________________
void ApplyListOfGatesOnAPureState(string s,MPS& psi,const SpinHalfSystem& C) {
  // The string s describes a list of operators/gates
  // format: op0_name q0a (q0b), op1_name q1a (q1b), ... 
  // |psi0> is replaced by |psi>=g0*g1*...*gN*|psi0>
  vector<string> ops=split(s,',');
  for (auto & op : ops) {
    vector<string> st=split(op,' ');
    string op_name;
    int i=-1,j=-1;
    int n=st.size();
    switch(n) {
      case 2:
        op_name=st[0];
        i=stoi(st[1]); if (i<1 || i>C.N) cout2<<"Error in SpinHalfSystem::ConstructProjectorFromGates: qubit index i="<<i<<" is out of range.\n",exit(0);
        if (op_name=="X" || op_name=="x") psi.ref(i)*=2*C.sites.op("Sx",i),psi.ref(i).noPrime();
        else if  (op_name=="Y" || op_name=="y") psi.ref(i)*=2*C.sites.op("Sy",i),psi.ref(i).noPrime();
        else if  (op_name=="Z" || op_name=="z") psi.ref(i)*=2*C.sites.op("Sz",i),psi.ref(i).noPrime();
        else if  (op_name=="H" || op_name=="h") psi.ref(i)*=Hadamard(C.sites,i),psi.ref(i).noPrime();
        else cout2<<"Error in SpinHalfSystem::ConstructProjectorFromGates: unknown 1-qubit operator "<<op_name<<".\n",exit(0);
        break;
      case 3:
        op_name=st[0];
        i=stoi(st[1]);
        j=stoi(st[2]);
        if (op_name=="CX" || op_name=="cx") ApplyControlledXYZGateONPureState(psi,C.sites,i,j,"Sx", Args("Cutoff",0));
        else if  (op_name=="CNOT" || op_name=="cnot") ApplyControlledXYZGateONPureState(psi,C.sites,i,j,"Sx", Args("Cutoff",0));
        else if  (op_name=="CY" || op_name=="cy") ApplyControlledXYZGateONPureState(psi,C.sites,i,j,"Sy", Args("Cutoff",0));
        else if  (op_name=="CZ" || op_name=="cz") ApplyControlledXYZGateONPureState(psi,C.sites,i,j,"Sz", Args("Cutoff",0));
        else cout2<<"Error in SpinHalfSystem::ConstructProjectorFromGates: unknown 2-qubit gate "<<op_name<<".\n",exit(0);
      break;
      default:
        cout2<<"Error in SpinHalfSystem::ConstructProjectorFromGates: expecting an operator name followed by 1 or 2 qubit number but got "<<op<<".\n",exit(0);
    }
  }
}

void StringToOperatorsList(string s, vector<string> &ops, vector<int> &qubits) {
  // The string s describes a list of operators
  // format: op0_name q0a, op1_name q1a, ...
  vector<string> l_ops=split(s,',');
  for (auto & op : l_ops) {
    vector<string> st=split(op,' ');
    string op_name;
    int i=-1;
    int n=st.size();
    switch(n) {
      case 2:
      {
        op_name=st[0];
        i=stoi(st[1]);
        if (i<1)  // || i>C.N)
            cout2 << "Error in StringToOperatorsList: qubit index i="<<i<<" is out of range.\n",exit(0);
        char op_lower = char(tolower(op_name[0]));
        if (op_lower=='x' || op_lower=='y' || op_lower=='z' || op_lower=='u' || op_lower=='d')
        {
            string s_op_lower = string("S");
            s_op_lower += op_lower;
            ops.push_back(s_op_lower);
            qubits.push_back(i);
        }
        else
            cout2 << "Error in StringToOperatorsList: unknown 1-qubit operator "<<op_name<<".\n",exit(0);
        break;
      }
      default:
        cout2 << "Error in StringToOperatorsList: expecting an operator name followed by 1 qubit number but got "<<op<<".\n",exit(0);
    }
  }
}