// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "grid3dconnectivity.h"
#include <iostream>

namespace dsl {

Grid3dConnectivity::Grid3dConnectivity(const Grid3d& grid)
    : LineConnectivity<Eigen::Vector3d, double >(grid) {
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        if (i == 0 && j == 0 && k == 0)
          continue;
        lines.push_back(Eigen::Vector3d(i, j, k));
      }
    }
  }
  for (auto& x : lines) {
    x = x.cwiseProduct(grid.cs);
    costs.push_back(x.norm());
  }
}

}
