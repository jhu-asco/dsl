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
#include "carcost.h"
#include "map.h"

namespace dsl {

using SE2Path = std::vector<Eigen::Matrix3d>;
using Vector1d = Eigen::Matrix<double, 1, 1>;

/**
 * Simple car connectivity using primitives. It enables the generation of arcs
 * from a given SE(2) state. These are then used to connect two cells in SE(2)
 * Each connection corresponds to an SE2Path, i.e. a vector of SE(2) poses
 *
 * Author: Marin Kobilarov
 */
class CarConnectivity : public GridConnectivity< Eigen::Vector3d, Eigen::Matrix3d, SE2Path > {
public:

  /**
   * Initialize cargrid connectivity with primitives corresponding to
   * motions with constant body fixed velocities for a fixed duration dt
   * @param grid The car grid
   * @param vs body-fixed velocities (each v=(vw,vx,vy))
   * @param dt time duration
   */
  CarConnectivity(const CarGrid& grid,
                  const std::vector<Eigen::Vector3d>& vs,
                  double dt = .25);
  
  
  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param onlyfwd If true, then only +ve forward velocity is used
   * @param wseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param tphimax Tan(max steering angle). Default value is 0.577=tan(M_PI/6)
   * @param vseg It decides the discretization of max forward velocity when
   * making the motion primitives
   * @param vxmax maximum forward velocity.
   */
  CarConnectivity(const CarGrid& grid,
                  bool onlyfwd = false,
                  int wseg = 4,
                  double tphimax = 0.577,
                  int vseg = 2,
                  double vxmax = 5);

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
                     int wseg = 4,
                     int vseg = 2);

  bool operator()(const SE2Cell& from,
                  std::vector< std::tuple<SE2Cell*, SE2Path, double> >& paths,
                  bool fwd = true) const override;

  bool Free(const Eigen::Matrix3d &g) const override { return true; }

  /**
   * Generate a path (a sequence of points) from initialize state in SE(2)
   * following
   * body-fixed velocity v=(vx,vy,w) for a unit time
   * @param path resulting path
   * @param g0 starting pose, matrix in SE(2)
   * @param v body fixed velocity (vx,vy,w)
   * @return true on success, false if obststructed by obstacle
   */
  bool Flow(std::tuple< SE2Cell*, SE2Path, double>& pathTuple,
            const Eigen::Matrix3d& g0,
            const Eigen::Vector3d& v) const;


  const CarGrid& grid; ///< the grid

  std::vector< Eigen::Vector3d >
      vs; ///< primitives defined using motions with constant
  /// body-fixed velocities (w,vx,vy)

  double dt = .5; ///< how long are the primitives

  CarCost cost;
};
}

#endif
