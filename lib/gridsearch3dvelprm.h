// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDSEARCH3DVELPRM_H
#define DSL_GRIDSEARCH3DVELPRM_H

#include <string.h>
#include "search.h"
#include "gridsearch3d.h"
#include "gridsearch3dvel.h"
#include "gridcost3dprm.h"
#include <iostream>

/**
 *  Grid based D* Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base Search class
 *
 *  The cell transition costs are hardcoded as the euclidean distance
 *  b/n the centers of the cells (i.e. each cell has 8 neighbors:
 *  the neighbors at N,S,W,E have transition cost of 1, and the 
 *  neighbors at NE, SE, NW, SW have transition costs of sqrt(2)
 *  These transition costs are added to the maximum of the values of
 *  two neighboring cells to compute the cost of their connecting edge.
 *  Finally that cost is multiplied by the variable scale.
 *
 *
 *  The planner is used as follows:
 *
 *  1) GridSearch(map, width, height)
 *  2) SetStart(x, y)
 *  3) SetGoal(x, y)
 *  4) Plan(path) and OptPath(path, optPath) -- to plan a path and optimize it
 *  5) follow path until map changes are observed
 *  6) for each change: SetCost(x, y, cost) or change the whole map: SetMap(map)
 *  7) SetStart(x,y) -- change the start to the current robot position
 *  8) goto 4
 * 
 *
 *  Author: Marin Kobilarov (c) 2004 mkobilar(at)robotics.usc.edu
 */

#define DSL3D_OCCUPIED 999999999

namespace dsl {


  class GridSearch3DVelPRM : public GridSearch3DVel
  {
  public:
    GridSearch3DVelPRM(int length, int width, int height, int numYaws, int numPitches, double v, double vmin[], double vmax[], double amin[], double amax[], const double *map = 0, double scale = 1.0) : GridSearch3DVel(length, width, height, numYaws, numPitches, map, scale), v(v)
    {
      cost.xub[0] = xub[0] = length;
      cost.xub[1] = xub[1] = width;
      cost.xub[2] = xub[2] = height;
      for(int i = 0; i < 3; i++)
      {
        cost.xlb[i] = xlb[i] = -5;
        cost.xlb[i+3] = xlb[i+3] = vmin[i];
        cost.xub[i+3] = xub[i+3] = vmax[i];
        cost.alb[i] = alb[i] = amin[i];
        cost.aub[i] = aub[i] = amax[i];
      }
    };
    
  protected:
    virtual void GetTrajectory(const Vertex &from, const Vertex &to, GridPath3DPlusTime &path) const;
 
    Graph graph;
    GridCost3DPRM cost; 
    Vertex **vertexMap;                ///< vertex grid array of size width*height
    double xlb[6];
    double xub[6];
    double alb[3];
    double aub[3];
    double v;
  };
}


#endif
