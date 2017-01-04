// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCOST_H
#define DSL_CARCOST_H

#include <memory>
#include <Eigen/Dense>
#include "gridcost.h"
#include "cargrid.h"
#include <memory.h>

namespace dsl {
/**
 * Car cost interface. Provides two cost interfaces: Heuristic and Real. Real returns the true cost of
 * moving from a cell to another cell based on some distance metric between 2 SE2Cell, also while checking
 * that the path along the twist is not blocked.
 *
 * Cost metrics: two cost metrics are implemented and is chosen via the use of corresponding constructor
 *   distance_metric =  ||pa-pb|| + ac*dist(aa, ab)), where pa,pb are the positions and aa,ab are the angles
 *   twistnorm_metric = ||weighted_twist||, where weighted_twist = w.diagonal*twist and wt is scaling weight
 */
  class CarCost : public GridCost< SE2Cell::PointType, SE2Cell::DataType > {
public:


  /**
   * Initialize the cost(distance_metric)
   * @param ac angular cost coefficient
   */
  CarCost(const CarGrid& grid, double ac, double eps = 1e-6);

  /**
   * Initialize the cost(twistnorm_metric)
   * @param wt weight for weighted norm of the twist between two SE2 cells
   */
  CarCost(const CarGrid& grid, const Vector3d& wt = Vector3d(0.1,1,2), double eps = 1e-6);

  /**
   * True cost(given a metric) of moving from cell "a" to a cell "b" along a twist element that connects them.
   * It also takes into account the occupancy structure of grid cells along the path and returns
   * numeric_limits<double>::quiet_NaN() cost if the path is blocked. It ban be checked with isnan() function
   *
   * @param a Start Cell
   * @param b End Cell
   * @return true cost of moving from a to b. If no path, returns numeric_limits<double>::quiet_NaN();
   */
  double Real(const SE2Cell& a, const SE2Cell& b) const;

  /**
   * Estimated cost(given a metric) of moving from cell "a" to a cell "b" along a twist element that connects them.
   * Unlike the cost returned by Real method, Heur doesn't check if path is blocked. The Heur cost is always an
   * underestimator of the Real cost.
   *
   * @param a Start cell
   * @param b End cell
   * @return estimated cost of moving from a to b.
   */
  double Heur(const SE2Cell& a, const SE2Cell& b) const;

  const CarGrid& grid_; ///< reference to the grid structure. Enables CarCost to give cost from cell a to b.

  bool use_twistnorm_metric_; ///< if true uses twistnorm_metric else distance_metric

  Vector3d wt_; ///< weight for weighted norm of the twist between two SE2 cells used by twistnorm_metric
  double ac_ = 1; ///< angular cost coefficient used by distance_metric(see explanation above),
  double eps_ = 1e-6; ///< for Heur cost a factor of (1-eps) is multiplied to final cost to ensure admissability
};
}

#endif
