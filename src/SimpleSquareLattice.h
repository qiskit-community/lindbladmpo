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

#ifndef _SIMPLESQUARELATTICE_
#define _SIMPLESQUARELATTICE_
//____________________________________________________________
//Super basic lattice class
//_____________________________________________________
class Lattice2d
{
public:
  int N;
  vector<int> X, Y; // Real space coordinates
  vector<int> I, J; //list of bonds
  map<pair<int, int>, int> index;
  bool predefined;
  bool predefined_chain;

  Lattice2d(){};
  //Constructor of some arbitrary lattice using two list of sites
  Lattice2d(int N_, vector<long> A, vector<long> B) : N(N_), predefined(false), predefined_chain(false)
  {
    cout << "A=" << A << "\nB=" << B << endl;
    if (N < 2)
      cerr << "Error: the number of sites N = " << N_ << ", must be >= 2.\n", exit(1);
    const unsigned int n = A.size();
    if (n != B.size())
      cerr << "Error: the two lists of sites have size " << A.size() << " and " << B.size() << " but they should be identical.\n", exit(1);
    //Loop over the bonds
    for (unsigned int i = 0; i < n; i++)
    {
      if (A[i] < 1 || A[i] > N)
        cerr << "Error, site number " << A[i] << " is not in [1,...,N=" << N << "].\n", exit(1);
      if (B[i] < 1 || B[i] > N)
        cerr << "Error, site number " << B[i] << " is not in [1,...,N=" << N << "].\n", exit(1);
      I.push_back(A[i]);
      J.push_back(B[i]);
    }
    cout << "User-defined lattice with " << N << " qbit(s) and " << n << " bonds.\n";
    cout << "List of bonds:\n";
    if (I.size() == 0)
    {
      cout << " [empty]" << endl;
    }
    else
    {
      cout << I << endl;
      cout << J << endl;
    }
  }
  // Creation of a Lx * Ly square lattice with or without periodic boundary conditions in the x and y direction
  Lattice2d(int Lx, int Ly, bool x_periodic = false, bool y_periodic = false) : N(Lx * Ly), predefined(true)
  {
    if (N < 2)
      cerr << "Error: the number of sites N = " << N << ", must be >= 2.\n", exit(1);
    predefined_chain = (Ly == 1) ? (true) : (false);
    bool up;
    int n = 0;
    for (int x = 0; x < Lx; x++)
    {
      if (x % 2 == 0)
        up = true;
      else
        up = false;
      if (up)
        for (int y = 0; y < Ly; y++)
        {
          X.push_back(x), Y.push_back(y);
          pair<int, int> p(x, y);
          index[p] = ++n;
        }
      if (!up)
        for (int y = Ly - 1; y >= 0; y--)
        {
          X.push_back(x), Y.push_back(y);
          pair<int, int> p(x, y);
          index[p] = ++n;
        }
    }
    for (auto it = index.cbegin(); it != index.cend(); ++it)
      std::cout << "x=" << it->first.first << " y=" << it->first.second
                << " #" << it->second << "\n";

    for (int x = 0; x < Lx; x++)
      for (int y = 0; y < Ly; y++)
      {
        pair<int, int> p(x, y);
        int ind = index[p];
        // "Vertical" bonds (parallel to the y axis)
        if (Ly > 1 && y + 1 < Ly)
        {
          I.push_back(ind);
          J.push_back(index[pair<int, int>(x, y + 1)]);
        }
        if (Ly > 2 && y_periodic == true && y == Ly - 1)
        {
          I.push_back(ind);
          J.push_back(index[pair<int, int>(x, 0)]);
        }
        // "Horizontal" bonds (parallel to the y axis)
        if (Lx > 1 && x + 1 < Lx)
        {
          I.push_back(ind);
          J.push_back(index[pair<int, int>(x + 1, y)]);
        }
        if (Lx > 2 && x_periodic == true && x == Lx - 1)
        {
          I.push_back(ind);
          J.push_back(index[pair<int, int>(0, y)]);
        }
      }
    cout << "Square lattice of size l_x=" << Lx << "*l_y=" << Ly << "=" << N << " qubits.\n";
    cout << "List of bonds:\n";
    if (I.size() == 0)
    {
      cout << " [empty]" << endl;
    }
    else
    {
      cout << I << endl;
      cout << J << endl;
    }
  }
};
//____________________________________________________________
#endif