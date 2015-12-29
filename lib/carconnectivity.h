// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID2DCONNECTIVITY_H
#define DSL_GRID2DCONNECTIVITY_H

#include "gridconnectivity.h"
#include "cargrid.h"
#include <vector>

namespace dsl {
  
  typedef GridPath<3, Matrix3d> SE2Path;

  /**
   * Simple car connectivity using primitives. It enables the generation of arcs
   * from a given SE(2) state. These are then used to connect two cells in SE(2)
   *
   * Author: Marin Kobilarov 
   */
  class CarConnectivity : public GridConnectivity<3, Matrix3d> {
  public:
    /**
     * Initialize a grid
     * @param grid the grid
     */
    CarConnectivity(const CarGrid &grid);

    /**
     * Use a set of primitive motions, i.e. arcs with body fixed forward velocities (-v,v) and 
     * angular velocity (-w,0,w), lasting time duration dt. There are 6 such combinations
     * @param vx forward velocity
     * @param w angular velocity
     * @param dt time duration
     */
    bool SetPrimitives(double vx, double w, double dt);
    
    bool operator()(const SE2Cell& from, 
                    std::vector<SE2Path>& paths,
                    bool fwd = true) const;

    /**
     * Generate a path (a sequence of points) from initialize state in SE(2) following
     * body-fixed velocity v=(vx,vy,w) for a unit time
     * @param path resulting path
     * @param g0 starting pose, matrix in SE(2)
     * @param v body fixed velocity (vx,vy,w)
     * @return true on success, false if obststructed by obstacle
     */
    bool Flow(SE2Path &path, const Eigen::Matrix3d &g0, const Eigen::Vector3d &v) const;
    
    const CarGrid &grid;    ///< the grid

    double vx;  ///< maximum forward velocity
    double w;   ///< maximum angular velocity
    double dt;  ///< how long are the primitives
    
    std::vector<Vector3d> vs;  ///< primitives defined using motions with constant body-fixed velocities (vx,vy,w)
  };
}

#endif
