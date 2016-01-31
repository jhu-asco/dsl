// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCONNECTIVITY_H
#define DSL_CARCONNECTIVITY_H

#include "gridconnectivity.h"
#include "cargrid.h"
#include <vector>

namespace dsl {

typedef GridPath< 3, Matrix3d > SE2Path;

/**
 * Simple car connectivity using primitives. It enables the generation of arcs
 * from a given SE(2) state. These are then used to connect two cells in SE(2)
 *
 * Author: Marin Kobilarov
 */
class CarConnectivity : public GridConnectivity< 3, Matrix3d > {
public:
  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param bp A scaling factor for penalizing going backward. Cost going
   * backward= bp*cost going forward
   * @param onlyfwd If true, then only +ve forward velocity is used
   * @param wseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param tphimax Tan(max steering angle). Default value is 0.577=tan(M_PI/6)
   * @param vseg It decides the discretization of max forward velocity when
   * making the motion primitives
   * @param vxmax maximum forward velocity.
   */
  CarConnectivity(const CarGrid& grid,
                  double bp = 1.0,
                  bool onlyfwd = false,
                  int wseg = 2,
                  double tphimax = 0.577,
                  int vseg = 1,
                  double vxmax = 1);

  /**
   * Use a set of primitive motions, i.e. arcs with body fixed forward
   * velocities (-v,v) and
   * angular velocity (-w,0,w), lasting time duration dt. By default there are 6
   * such combinations.
   * In addition we made it configurable so that more primitves with angular
   * velocities in between(-w,w)
   * can be added in between by choosing the wseg parameter >2. Also only +ve
   * forward velocity can be chosen
   * @param vx forward velocity
   * @param w angular velocity
   * @param dt time duration
   * @param onlyfwd If true, then primitives with +ve forward velocity are made
   * @param wseg The discretization of max angular velocity for making the
   * primitives
   * @param vseg The discretization of max forward velocity for making the
   * primitives
   * @return
   */
  bool SetPrimitives(double vx,
                     double w,
                     double dt,
                     double onlyfwd = false,
                     int wseg = 2,
                     int vseg = 1);

  bool operator()(const SE2Cell& from,
                  std::vector< SE2Path >& paths,
                  bool fwd = true) const;

  /**
   * Generate a path (a sequence of points) from initialize state in SE(2)
   * following
   * body-fixed velocity v=(vx,vy,w) for a unit time
   * @param path resulting path
   * @param g0 starting pose, matrix in SE(2)
   * @param v body fixed velocity (vx,vy,w)
   * @return true on success, false if obststructed by obstacle
   */
  bool Flow(SE2Path& path,
            const Eigen::Matrix3d& g0,
            const Eigen::Vector3d& v) const;

  const CarGrid& grid; ///< the grid

  double vx; ///< maximum forward velocity
  double w;  ///< maximum angular velocity
  double dt; ///< how long are the primitives
  double bp; ///< A scaling factor for penalizing going backward. Cost going
  /// backward= bp*cost going forward

  std::vector< Vector3d >
      vs; ///< primitives defined using motions with constant
  /// body-fixed velocities (w,vx,vy)
};
}

#endif
