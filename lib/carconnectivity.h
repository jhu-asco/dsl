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
  
  /**
   * Simple car connectivity using primitives. It enables the generation of arcs
   * from a given SE(2) state. These are then used to connect two cells in SE(2)
   *
   * Author: Marin Kobilarov 
   */
  class CarConnectivity : public GridConnectivity<3> {
  public:
    /**
     * Initialize a grid
     * @param grid the grid
     */
    CarConnectivity(const CarGrid &grid);
    
    bool operator()(const Cell<3>& from, 
                    std::vector<GridPath<3> >& paths,
                    bool fwd = true) const;

    /**
     * Generate a path (a sequence of points) from initialize state in SE(2) following
     * body-fixed velocity v=(vx,vy,w) for a unit time
     */
    bool Flow(GridPath<3> &path, const Eigen::Matrix3d &g0, const Eigen::Vector3d &v) const;
    
    const CarGrid &grid;    ///< the grid

    double v;  ///< maximum forward velocity
    double w;  ///< maximum angular velocity
    double dt; ///< how long are the primitives
    
    std::vector<Vector3d> vs;  ///< primitives defined using motions with constant body-fixed velocities (vx,vy,w)
  };
}

#endif
