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

  using SE2Cell =  Cell< Eigen::Vector3d, Eigen::Matrix3d >;
  using SE2CellPtr = shared_ptr<SE2Cell>;

// geometry of the car
struct CarGeom {
  CarGeom (double l = 3, double b = 1.5, double ox = 0.75, double oy = 0);
  double l;  ///< dim of the rectangle bounding the car along the x direction (from back to front)
  double b;  ///< dim of the rectangle bounding the car along the y direction (from right to left)
  double ox; ///< x position of the center of the origin of the car wrt to the center of bounding rectangle
  double oy; ///< y position of the center of the origin of the car wrt to the center of bounding rectangle

  /**
   * Fill the rectange representing the car with regularly spaced points
   * @param cs
   * @param points
   */
  void Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const;
};


/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell also
 *stores a 3x3 SE(2) matrix describing the pose
 * of the car.
 *
 * Authors: Marin Kobilarov, Subhransu Mishra
 */
class CarGrid : public Grid< Eigen::Vector3d, Eigen::Matrix3d > {
public:
  CarGrid(const Map<bool, 3> &cmap,
          const Eigen::Vector3d& cs);
  
  const Map<bool, 3>& cmap; ///< configuration-space map
  
  static void MakeMap(const Map<bool, 2> &map, Map<bool, 3> &cmap);
  
  static void MakeMap(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap);

  static void DilateMap(const CarGeom& geom, double theta,
                          double sx, double sy, int gx, int gy,
                          const vector<bool>& data, vector<bool>& data_dil);

};
}

#endif
