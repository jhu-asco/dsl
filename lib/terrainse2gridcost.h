// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_TERRAINSE2GRIDCOST_H
#define DSL_TERRAINSE2GRIDCOST_H

#include <Eigen/Dense>
#include "gridcost.h"
#include <assert.h>
#include <terrainse2grid.h>

namespace dsl {

/**
 * Cost interface for a 3D grid over SE2 pose, parametrized by yaw, x and y, where each
 * grid cells hold the terrain properties like height and traversibility.Provides two cost
 * interfaces: Heuristic and Real. Real returns the true cost of moving from a cell to
 * another cell based on some distance metric between 2 SE2Cell, also while checking that
 * the path along the twist is not blocked.
 *
 */
class TerrainSE2GridCost : public GridCost< TerrainCell::PointType, TerrainCell::DataType > {
public:
  /**
   * Initialize the cost(twistnorm_metric)
   * @param wt weight for weighted norm of the twist between two SE2 cells
   */
  TerrainSE2GridCost(const TerrainSE2Grid& grid, const Vector3d& wt = Vector3d(0.1,1,2), double eps = 1e-6);

  /**
   * True cost(given a metric) of moving from cell "a" to a cell "b" along a twist element that connects them.
   * It also takes into account the occupancy structure of grid cells along the path and returns
   * numeric_limits<double>::quiet_NaN() cost if the path is blocked. It ban be checked with isnan() function
   *
   * @param a Start Cell
   * @param b End Cell
   * @return true cost of moving from a to b. If no path, returns numeric_limits<double>::max();
   */
  double Real(const TerrainCell& a, const TerrainCell& b) const;

  /**
   * Estimated cost(given a metric) of moving from cell "a" to a cell "b" along a twist element that connects them.
   * Unlike the cost returned by Real method, Heur doesn't check if path is blocked. The Heur cost is always an
   * underestimator of the Real cost.
   *
   * @param a Start cell
   * @param b End cell
   * @return estimated cost of moving from a to b.
   */
  double Heur(const TerrainCell& a, const TerrainCell& b) const;

  const TerrainSE2Grid& grid; ///< reference to the grid structure. Enables CarCost to give cost from cell a to b.
  Vector3d wt; ///< weight for weighted norm of the twist between two SE2 cells used by twistnorm_metric
  double eps = 1e-6; ///< for Heur cost a factor of (1-eps) is multiplied to final cost to ensure admissability
};
}

#endif // DSL_TERRAINSE2GRIDCOST_H
