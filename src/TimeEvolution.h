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

#ifndef _TIMEEVOLUTION_
#define _TIMEEVOLUTION_

#include "itensor/all.h"
#include <string>

using namespace itensor;
using namespace std;

//____________________________________________________________________
class TimeEvolver
{
  public:
    // The MPO below are associated to a given time step tau and an given Lindbladian
    // They are used to implement the time evolution.
    int order; // 2,3 or 4
    MPO expL1, expL2, expL3, expL4, expL5, expL6, expL7;

    // Object to store possible options related to the time evolution
    Args argsApplyMPOtoRho;

    // The 'init' below has be be called once, so that the expL1...expL7 above are constructed
    // if the Lindbladian and/or the time step changes, then init has to be called again
    void init(double tau, const AutoMPO &auto_L, Args args, int ord = 4);
    // Actual time evolution (1 'small' time step tau [value defined])
    void evolve(MPS &rho) const;
};
//____________________________________________________________________
#endif