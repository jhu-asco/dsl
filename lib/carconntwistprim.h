// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu> and Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCONNTWISTPRIM_H
#define DSL_CARCONNTWISTPRIM_H

#include "gridconnectivity.h"
#include "cargrid.h"
#include <vector>
#include "carcost.h"
#include "map.h"
#include "carprimcfg.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace dsl {

using SE2Prim = Eigen::Vector3d;
using Vector1d = Eigen::Matrix<double, 1, 1>;

/**
 * Simple car connectivity using primitives. It enables the generation of arcs
 * from a given SE(2) state. These are then used to connect two cells in SE(2)
 * Each connection corresponds to an SE2Prim, i.e. a vector of SE(2) poses
 *
 * Author: Marin Kobilarov
 */
class CarConnTwistPrim: public GridConnectivity< Eigen::Vector3d, Eigen::Matrix3d, SE2Prim > {
public:

  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param cfg The configuration for generating primitives
   */
  CarConnTwistPrim(const CarGrid& grid, CarPrimitiveCfg& cfg);

  /**
   * Initialize cargrid connectivity with primitives corresponding to
   * motions with constant body fixed velocities for a fixed duration dt
   * @param grid The car grid
   * @param vs body-fixed velocities (each v=(vw,vx,vy))
   * @param dt time duration
   */
  CarConnTwistPrim(const CarGrid& grid,
                  const std::vector<Eigen::Vector3d>& vs,
                  double dt = .25);
  
  
  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param dt time duration
   * @param vx forward velocity.
   * @param kmax maximum curvature k = Tan(max steering angle)/axle_length; Default value is 0.577=tan(M_PI/6)/1
   * @param kseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param onlyfwd If true, then only +ve forward velocity is used
   */
  CarConnTwistPrim(const CarGrid& grid,
                  double dt = .25,
                  double vx = 5,
                  double kmax = 0.577,
                  int kseg = 4,
                  bool onlyfwd = false);
  
  /**
   * Generates primitives based on the configuration
   * @param dt time duration
   * @param cfg the configuration for generation of primitives
   * @return
   */
  bool SetPrimitives(double dt, CarPrimitiveCfg& cfg);

  /**
   * Use a set of primitive motions, i.e. arcs with body fixed forward
   * velocities (-v,v) and
   * angular velocity (-w,0,w), lasting time duration dt. By default there are 6
   * such combinations.
   * In addition we made it configurable so that more primitves with angular
   * velocities in between(-w,w)
   * can be added in between by choosing the kseg parameter >2. Also only +ve
   * forward velocity can be chosen
   * @param dt time duration
   * @param vx forward velocity.
   * @param kmax maximum curvature k = Tan(max steering angle)/axle_length; Default value is 0.577=tan(M_PI/6)/1
   * @param kseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param onlyfwd If true, then only +ve forward velocity is used
   * @return true on success
   */
  bool SetPrimitives(double dt,
                     double vx,
                     double kmax,
                     int kseg = 4,
                     bool onlyfwd = false);
  
  bool operator()(const SE2Cell& from,
                  std::vector< std::tuple<SE2Cell*, SE2Prim, double> >& paths,
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
  bool Flow(std::tuple< SE2Cell*, SE2Prim, double>& pathTuple,
            const Eigen::Matrix3d& g0,
            const Eigen::Vector3d& v) const;


  const CarGrid& grid; ///< the grid

  std::vector< Eigen::Vector3d > vs; ///< primitives defined using motions with constant body-fixed velocities (w,vx,vy)
  std::vector< std::vector<Eigen::Vector3d > > vss; ///< A set of primitives for each angle

  double dt = .5; ///< how long are the primitives

  CarCost cost;
};
}

#endif //DSL_CARCONNTWISTPRIM_H
