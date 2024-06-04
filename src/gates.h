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

#ifndef _GATES_
#define _GATES_

#include "Pauli.h"
#include "io_utils.h"
#include <string>

// Apply the X gate (qubit i) on a mixed state rho
void ApplyXGate(MPS &, const Pauli &, int);

// Apply the Y gate (qubit i) on a mixed state rho
void ApplyYGate(MPS &, const Pauli &, int);

// Apply the X gate (qubit i) on a mixed state rho
void ApplyZGate(MPS &r, const Pauli &, int);

// Apply the SqrtX gate (qubit i) on a mixed state rho
void ApplySqrtXGate(MPS &, const Pauli &, int);

void ApplyHGate(MPS &, const Pauli &, int);

void ApplyProjUp(MPS &, const Pauli &, int);

void ApplyProjDn(MPS &, const Pauli &, int);

void ApplyControlledXYZGate(MPS &, const Pauli &, int, int, string, Args = Args("Cutoff", 0));

// Apply the CNOT gate on some mixed state rho, at sites (control,target)
void ApplyCNOTGate(MPS &, const Pauli &, int, int, Args = Args("Cutoff", 0));

// Apply the controlled-Z gate on some mixed state rho, at sites (i,j)
void ApplyControlledZGate(MPS &, const Pauli &, int, int, Args arg = Args("Cutoff", 0));

// From an ordered list of gates [g1,g2, ..., gN] specified in a string s, and from a 'initial' pure state |psi0>,
// construct |psi>=g1*g2*...*gN |psi0>
void ApplyListOfGatesOnAPureState(string, MPS &, const SpinHalfSystem &);

void StringToOperatorsList(string, vector<string> &, vector<int> &);

#endif