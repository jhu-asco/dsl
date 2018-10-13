#include <iostream>
#include "grid2dconnectivity.h"

namespace dsl {

Grid2dConnectivity::Grid2dConnectivity(
  const Grid2d& grid)
    : LineConnectivity<Eigen::Vector2d, double >(grid) {
  lines.push_back(Eigen::Vector2d(-1, -1));
  lines.push_back(Eigen::Vector2d(0, -1));
  lines.push_back(Eigen::Vector2d(1, -1));
  lines.push_back(Eigen::Vector2d(-1, 0));
  lines.push_back(Eigen::Vector2d(1, 0));
  lines.push_back(Eigen::Vector2d(-1, 1));
  lines.push_back(Eigen::Vector2d(0, 1));
  lines.push_back(Eigen::Vector2d(1, 1));

  for (auto&& x : lines) {
    x = x.cwiseProduct(grid.cs); // scale by grid cell size
    costs.push_back(x.norm());
  }
}

}
