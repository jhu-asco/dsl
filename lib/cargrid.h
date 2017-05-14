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
#include "map.h"

namespace dsl {

  typedef Cell< Eigen::Vector3d, Eigen::Matrix3d > SE2Cell;

// geometry of the car
struct CarGeom {
  CarGeom (double l = 0, double b = 0, double ox = 0, double oy = 0) : l(l), b(b), ox(ox), oy(oy) {}
  double l;  ///< dim of the rectangle bounding the car along the x direction
  double b;  ///< dim of the rectangle bounding the car along the y direction

  //  Eigen::Vector2d o;
  
  double ox; ///< x position of the center of the origin of the car wrt to the
  /// center of bounding rectangle
  double oy; ///< y position of the center of the origin of the car wrt to the

  void Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const {
    points.clear();
    for (double x = 0; x <= l; x += cs[0])
      for (double y = 0; y <= b; y += cs[1])
        points.push_back(Eigen::Vector2d(x + ox, y + oy));
  }

  /// center of bounding rectangle
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
  
  static void MakeMap(const CarGeom& geom, const Map<bool, 2> &map, Map<bool, 3> &cmap);
  
  static void DilateMap(const CarGeom& geom, double theta,
                        double sx, double sy, int gx, int gy, 
                        const bool* data, bool* data_dil);


  void Slice(const Map<bool, 3> &cmap, double a, Map<bool, 2> &map) const;

  
    //  std::vector<Eigen::Vector2d> raster;
};
}

#endif
