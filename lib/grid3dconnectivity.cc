#include "grid3dconnectivity.h"
#include <iostream>

namespace dsl {

using Eigen::Vector3d;

Grid3dConnectivity::Grid3dConnectivity(const Grid< Vector3d, double >& grid)
    : LineConnectivity< Vector3d, double >(grid) {
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        if (i == 0 && j == 0 && k == 0)
          continue;
        lines.push_back(Vector3d(i, j, k));
      }
    }
  }
  for (auto& x : lines) {
    x = x.cwiseProduct(grid.cs());
    costs.push_back(x.norm());
  }
}

}
