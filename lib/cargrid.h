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
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>
namespace dsl {

typedef Cell<3, Matrix3d> SE2Cell;

class CarGrid : public Grid<3, Matrix3d> {
public:

  CarGrid(double l,double b, double ox, double oy,
          int width, int height, double* map, double sx, double sy, double sa,
          double costScale,double maxCost=1);

  CarGrid(int width, int height, double *map,
          double sx, double sy, double sa, double costScale,
          double maxCost = 1);

  double maxCost; ///< any cell cost above maxCost is considered obstacle and not added to the graph
  double l_;      ///< dim of the rectange bounding the car along the x direction
  double b_;      ///< dim of the rectange bounding the car along the y direction
  double ox_;     ///< x position of the center of the origin of the car wrt to the center of bounding rectange
  double oy_;     ///< y position of the center of the origin of the car wrt to the center of bounding rectange

private:
  void getDilatedMap(double* data_dil, double* data, double theta);

};
}

#endif
