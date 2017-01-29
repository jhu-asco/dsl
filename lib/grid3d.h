// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GRID3D_H_
#define DSL_LIB_GRID3D_H_

#include "grid.h"

namespace dsl {

using XyzCostCell = Cell<Eigen::Vector3d, double>;
using Grid3dBase = Grid<XyzCostCell::PointType, XyzCostCell::Ptr>;
/**
 * A 3d grid with coordinates (x,y,z)
 *
 * Author: Matt Sheckells
 */
  class Grid3d : public Grid3dBase {
public:
  /**
   * Initialize the grid using a 3d configuration-space map
   * @param length length
   * @param width width
   * @param height height
   * @param map the map
   * @param sx x-scale of each grid cell
   * @param sy y-scale of each grid cell
   * @param sz z-scale of each grid cell
   * @param costScale multiple the value of each grid cell by costScale
   * @param maxCost any cell cost above maxCost is considered obstacle and not
   * added to the graph
   */
  Grid3d(int length,
         int width,
         int height,
         const double* map,
         double sx,
         double sy,
         double sz,
         double costScale,
         double maxCost = 1);
};
}

#endif
