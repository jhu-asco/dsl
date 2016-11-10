// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCOST_H
#define DSL_CARCOST_H

#include <Eigen/Dense>
#include "gridcost.h"
#include <memory>

namespace dsl {
using std::shared_ptr;

  // a cell that stores an SE(2) transformation  matrix
  using SE2Cell =  Cell< Eigen::Vector3d, Eigen::Matrix3d >;
  using SE2CellPtr = shared_ptr<SE2Cell>;

/**
 * Car cost interface
 *
 */
  class CarCost : public GridCost< Eigen::Vector3d, Eigen::Matrix3d > {
public:
  /**
   * Initialize the cost
   * @param ac angular cost coefficient
   */
  CarCost(double ac = 1, double eps = 1e-6);

  double Heur(const SE2Cell& a, const SE2Cell& b) const;
  double Real(const SE2Cell& a, const SE2Cell& b) const;

  double ac = 1; ///< angular cost coefficient, default is 1 (the total cost is
  /// proportional to |pa-pb| + ac*dist(aa, ab) ), where pa,pb are the
  /// positions and aa,ab are the angles

  double eps = 1e-6; ///< heuristic cost is some quadratic norm on the coordinates multiplied by (1-eps)
};
}

#endif
