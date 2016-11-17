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

  void set(double l, double b, double ox, double oy, double sb = 0);

  /**
   * Fill the rectangle representing the car with regularly spaced points
   * @param cs
   * @param points
   */
  void Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const;

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
  
  /**
   * Takes a 2-dimensional occupancy grid and repeats that for all angles
   * @param omap  2-dimensional occupancy grid
   * @param cmap  3-dimensional occupancy grid
   * @return True if cmap could be modified correctly
   */
  static bool MakeMap(const Map<bool, 2> &omap, Map<bool, 3> &cmap);
  
  /**
   * Takes a 2-dim occupancy grid and uses car geometry to create a 3-dim occ grid for all angles
   * @param geom Geometry of the car
   * @param omap 2-dimensional occupancy grid
   * @param cmap 3-dimensional occupancy grid
   * @return True if cmap could be modified correctly
   */
  static bool MakeMap(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap);

  static void DilateMap(const CarGeom& geom, double theta,
                          double sx, double sy, int gx, int gy,
                          const vector<bool>& data, vector<bool>& data_dil);

};
}

#endif
