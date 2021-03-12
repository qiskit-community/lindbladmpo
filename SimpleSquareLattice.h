#ifndef _SIMPLESQUARELATTICE_
#define _SIMPLESQUARELATTICE_
//____________________________________________________________
//Super basic square lattice class, to keep track of the site positions, bonds, etc.
//_____________________________________________________
class Lattice2d
{
public:
  int N;
  vector<int> X, Y; // Real space coordinates
  vector<int> I, J; //list of bonds
  map<pair<int, int>, int> index;
  Lattice2d(int Lx, int Ly, bool x_periodic = false, bool y_periodic = false) : N(Lx * Ly)
  {
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
          J.push_back(index[pair<int, int>(x , y+1)]);
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
    cout << "Bonds:\n";
    cout << I << endl;
    cout << J << endl;
  }
};
//____________________________________________________________
#endif