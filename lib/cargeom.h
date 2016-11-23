// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARGEOM_H
#define DSL_CARGEOM_H

#include "grid.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

namespace dsl {

  using Vector5d = Eigen::Matrix<double,5,1>;

/**
 * Geometry of a car defined as a rectangle with some origin. The position of the origin is
 * defined relative to the center of rectangle. Safety buffer basically grows the rectange
 * outward to added safety. For example if you consider a car with length 3m,  width 1.5m,
 * origin at middle point of rear axle and  safety buffer of 40cm, then we have
 * l = 3.0, b = 1.5, ox = 0.75, oy = 0.0, sb = 0.4
 * TODO: Add constructors with a kernel image to indicate the shape of car
 */
class CarGeom {
public:
  CarGeom (double l = 0, double b = 0, double ox = 0, double oy = 0, double sb = 0);

  CarGeom (const Vector5d& lboxoysb);

  void set(double l, double b, double ox, double oy, double sb = 0);

  void set(const Vector5d& lboxoysb);
  /**
   * Fill the rectangle representing the car with regularly spaced points
   * @param cs
   * @param points
   */
  void Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const;

  /**
   * Get Corners of the car in a cyclic fashion
   * @param vertices
   * @param theta
   */
  void GetTrueCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const;

  /**
   * Get Corners of the car(with safety buffer) in a counter clockwise fashion starting
   * @param vertices
   * @param theta
   */
  void GetSafeCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const;

  double l() const;
  double b() const;
  double ox() const;
  double oy() const;
  double sb() const;
  double le() const;
  double be() const;

private:
  double l_;  ///< dim of the rectangle bounding the car along the x direction (from back to front)
  double b_;  ///< dim of the rectangle bounding the car along the y direction (from right to left)
  double ox_; ///< x position of the center of the origin of the car wrt to the center of bounding rectangle
  double oy_; ///< y position of the center of the origin of the car wrt to the center of bounding rectangle
  double sb_; ///< safety buffer. The distance to which the car rectange is expanded

  double le_; ///<Effective l because of safety buffer
  double be_; ///<Effective b because of safety buffer
};
}

#endif /* DSL_CARGEOM_H */
