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
