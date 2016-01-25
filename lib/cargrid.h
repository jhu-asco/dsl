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

typedef Cell< 3, Matrix3d > SE2Cell;

// geometry of the car
struct CarGeom {
  CarGeom (double l = 0, double b = 0, double ox = 0, double oy = 0) : l(l), b(b), ox(ox), oy(oy) {}
  double l;  ///< dim of the rectangle bounding the car along the x direction
  double b;  ///< dim of the rectangle bounding the car along the y direction
  double ox; ///< x position of the center of the origin of the car wrt to the
  /// center of bounding rectangle
  double oy; ///< y position of the center of the origin of the car wrt to the
  /// center of bounding rectangle
};

// 2d map of the environment
struct Map2d {
  Map2d(int width = 0, int height = 0, const double *data = 0) :
    width(width), height(height), data(data) {}
  int width;
  int height;
  const double* data;   
};

/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell also
 *stores a 3x3 SE(2) matrix describing the pose
 * of the car.
 *
 * Authors: Marin Kobilarov, Subhransu Mishra
 */
class CarGrid : public Grid< 3, Matrix3d > {
public:
  CarGrid(const CarGeom &geom,
          const Map2d &map,
          double sx,
          double sy,
          double sa,
          double costScale,
          double maxCost = 1);

  /**
   * Initialize the grid using a 2d configuration-space map
   * @param width width
   * @param height height
   * @param map the map
   * @param sx x-scale of each grid cell
   * @param sy y-scale of each grid cell
   * @param sa angle-scale of each grid cell
   * @param costScale multiple the value of each grid cell by costScale
   * @param maxCost any cell cost above maxCost is considered obstacle and not
   * added to the graph
   */
  CarGrid(const Map2d& map,
          double sx,
          double sy,
          double sa,
          double costScale,
          double maxCost = 1);

  double maxCost; ///< any cell cost above maxCost is considered obstacle and
  /// not added to the graph

  CarGeom geom;  ///< car geometry

private:
  void getDilatedMap(double* data_dil, const double* data, double theta);
};
}

#endif
