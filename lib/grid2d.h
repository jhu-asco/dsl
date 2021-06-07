// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GRID2D_H_
#define DSL_LIB_GRID2D_H_

#include "grid.h"

namespace dsl {

using XyCostCell = Cell<Eigen::Vector2d, double>;
using Grid2dBase = Grid<XyCostCell::PointType, XyCostCell>;

/**
 * A 2d grid with coordinates (x,y)
 *
 */
  class Grid2d : public Grid2dBase {
public:
  /**
   * Initialize the grid using a 2d configuration-space map
   * @param width width
   * @param height height
   * @param map the map
   * @param sx x-scale of each grid cell
   * @param sy y-scale of each grid cell
   * @param maxCost any cell cost above maxCost is considered obstacle and not
   * added to the graph
   */
  Grid2d(int width,
         int height,
         const double* map,
         double sx,
         double sy,
         double maxCost = 1);
};
}

#endif
