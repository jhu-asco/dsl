// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GRID3DCONNECTIVITY_H_
#define DSL_LIB_GRID3DCONNECTIVITY_H_

#include "lineconnectivity.h"

namespace dsl {

/**
 * Defines a simple connectivity between cells in a 3d grid.
 * The default implementation is the 26-cell Moore neighborhood connectivity.
 * The costs are the Euclidean distances b/n the cell centers.
 */
  class Grid3dConnectivity : public LineConnectivity< Eigen::Vector3d, double > {
public:
    Grid3dConnectivity(const Grid< Eigen::Vector3d, double >& grid);
};
}

#endif
