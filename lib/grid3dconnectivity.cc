#include "grid3dconnectivity.h"
#include <iostream>

using namespace dsl;
using namespace std;

Grid3dConnectivity::Grid3dConnectivity(const Grid<3> &grid) : LineConnectivity<3>(grid) {
  for(int i = -1; i <= 1; i++)
  {
    for(int j = -1; j <= 1; j++)
    {
      for(int k = -1; k <= 1; k++)
      {
        if(i == 0 && j == 0 && k == 0)
          continue;
        lines.push_back(Vector3d(i,j,k));                          
      }
    }
  }
  vector<Vectornd>::iterator it;
  for (it = lines.begin(); it != lines.end(); ++it) {
    Vectornd &x = *it;
    x = x.cwiseProduct(grid.cs);
    costs.push_back(x.norm());
  }
}
