// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "lineconnectivity.h"
#include "grid3d.h"

namespace dsl {

/**
 * Defines a simple connectivity between cells in a 3d grid.
 * The default implementation is the 26-cell Moore neighborhood connectivity.
 * The costs are the Euclidean distances b/n the cell centrs.
 */
  class Grid3dConnectivity : public LineConnectivity< Eigen::Vector3d, double > {
public:
    Grid3dConnectivity(const Grid3d& grid);
};
}
