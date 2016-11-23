// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARGRID_H
#define DSL_CARGRID_H

#include "grid.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

namespace dsl {
using std::shared_ptr;

using SE2Cell =  Cell< Eigen::Vector3d, Eigen::Matrix3d >; // a cell that stores an SE(2) transformation  matrix
using SE2CellPtr = shared_ptr<SE2Cell>;
using SE2CellGrid = Grid<SE2Cell::PointType, SE2Cell::DataType>;

/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell also
 *stores a 3x3 SE(2) matrix describing the pose
 * of the car.
 *
 * Authors: Marin Kobilarov, Subhransu Mishra
 */
class CarGrid : public Grid< SE2Cell::PointType, SE2Cell::DataType > {
public:

  CarGrid(const Map<bool, 3> &cmap, const Eigen::Vector3d& cs);
  const Map<bool, 3>& cmap; ///< configuration-space map
};
}

#endif
