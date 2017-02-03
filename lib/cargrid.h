// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_CARGRID_H_
#define DSL_LIB_CARGRID_H_

#include "grid.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

namespace dsl {
using SE2Cell =  Cell< Eigen::Vector3d, Eigen::Matrix3d >; // a cell that stores an SE(2) transformation  matrix
using SE2CellGrid = Grid<SE2Cell::PointType, SE2Cell>;
/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell also
 *stores a 3x3 SE(2) matrix describing the pose
 * of the car.
 *
 * Authors: Marin Kobilarov, Subhransu Mishra
 */
class CarGrid : public SE2CellGrid {
public:

  /**
   * Creates the lattice grid and allocates memory for cells that are free.
   * @param cmap configuration-space map(occupancy grid in 3-dim)
   * @param cs the cell size to use for the grid. Generally CarGrid.cs > cmap.cs
   */
  CarGrid(const Map<bool, 3> &cmap, const Eigen::Vector3d& cs);

};
}

#endif
