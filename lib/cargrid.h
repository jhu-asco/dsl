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

typedef Cell<3, Matrix3d> SE2Cell;

/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell also stores a 3x3 SE(2) matrix describing the pose
 * of the car.
 *
 * Authors: Marin Kobilarov, Subhransu Mishra
 */
class CarGrid : public Grid<3, Matrix3d> {
public:

  CarGrid(double l,double b, double ox, double oy,
          int width, int height, const double* map, double sx, double sy, double sa,
          double costScale,double maxCost=1);

  /**
   * Initialize the grid using a 2d configuration-space map
   * @param width width
   * @param height height
   * @param map the map
   * @param sx x-scale of each grid cell
   * @param sy y-scale of each grid cell
   * @param sa angle-scale of each grid cell
   * @param costScale multiple the value of each grid cell by costScale
   * @param maxCost any cell cost above maxCost is considered obstacle and not added to the graph
   */
  CarGrid(int width, int height, const double *map,
          double sx, double sy, double sa, double costScale,
          double maxCost = 1);

  double maxCost; ///< any cell cost above maxCost is considered obstacle and not added to the graph
  double l_;      ///< dim of the rectangle bounding the car along the x direction
  double b_;      ///< dim of the rectangle bounding the car along the y direction
  double ox_;     ///< x position of the center of the origin of the car wrt to the center of bounding rectangle
  double oy_;     ///< y position of the center of the origin of the car wrt to the center of bounding rectangle

private:

  void getDilatedMap(double* data_dil, const double* data, double theta);

};
}

#endif
