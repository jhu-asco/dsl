#include "grid2dconnectivity.h"
#include <iostream>

using namespace dsl;
using namespace std;

Grid2dConnectivity::Grid2dConnectivity(const Grid<2> &grid) : LineConnectivity<2>(grid) {
  
  lines.push_back(Vector2d(-1, -1));
  lines.push_back(Vector2d(0, -1));
  lines.push_back(Vector2d(1, -1));
  lines.push_back(Vector2d(-1, 0));
  lines.push_back(Vector2d(1, 0));
  lines.push_back(Vector2d(-1, 1));
  lines.push_back(Vector2d(0, 1));
  lines.push_back(Vector2d(1, 1));                          

  vector<Vectornd>::iterator it;
  for (it = lines.begin(); it != lines.end(); ++it) {
    Vectornd &x = *it;
    x = x.cwiseProduct(grid.cs);
    costs.push_back(x.norm());
  }
}
