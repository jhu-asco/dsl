#include "grid2dconnectivity.h"
#include <iostream>

namespace dsl {

using Eigen::Vector2d;

Grid2dConnectivity::Grid2dConnectivity(const Grid2dBase& grid)
    : Line2dConnectivity(grid) {
  lines.push_back(Vector2d(-1, -1));
  lines.push_back(Vector2d(0, -1));
  lines.push_back(Vector2d(1, -1));
  lines.push_back(Vector2d(-1, 0));
  lines.push_back(Vector2d(1, 0));
  lines.push_back(Vector2d(-1, 1));
  lines.push_back(Vector2d(0, 1));
  lines.push_back(Vector2d(1, 1));
  
  for (auto&& x : lines) {
    x = x.cwiseProduct(grid.cs()); // scale by grid cell size
    costs.push_back(x.norm());
  }
}

}
