// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GEOM_H_
#define DSL_LIB_GEOM_H_

#include <Eigen/Dense>
#include <vector>

namespace dsl {

/**
 * Interface defining the geometry of the robot in dimension
 * The geometry should encode the physical structure of the robot as well as the desired
 * clearance during planning.
 *
 * This interface defines method to access the geometry of the robot in different ways:
 * * set of corners
 * * 2D or 3D grid of occupancy
 *
 * Author: Subhransu Mishra -- Copyright (C) 2004
 */
template < class PointType>
class Geom{
  /**
   * Fill the rectangle representing the car with regularly spaced points
   * @param cs
   * @param points
   */
  virtual void Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const = 0;

  /**
   * Get Corners of the car in a cyclic fashion
   * @param vertices
   * @param theta
   */
  virtual void GetTrueCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const = 0;

  /**
   * Get Corners of the car(with safety buffer) in a counter clockwise fashion starting
   * @param vertices
   * @param theta
   */
  virtual void GetSafeCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const = 0;
};


} /* namespace dsl */

#endif /* DSL_LIB_GEOM_H_ */
