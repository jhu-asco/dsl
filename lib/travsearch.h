// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_TRAVSEARCH_H
#define DSL_TRAVSEARCH_H

#include "gridsearch.h"
#include "gridcost.h"

/**
 *  Grid based D* Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base Search class.  The cost function in particular
 *  uses the gradient of the vertex costs across an edge length to define the edge cost.
 *  The purpose of this is to generate a traversability map where the vertex costs
 *  are the height of a particular 2D position.  
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
 *  Author: Matt Sheckells (c) 2015 msheckells(at)jhu.edu
 */


namespace dsl {

  class TravSearch : public GridSearch
  {
  public:
    
    /**
     * The planner allocates all states in the beginning
     * a valid map of size width*height should be provided
     * scale is used to multiply each cell value and each neighbor transition value
     * before they are used for planning
     * @param width width
     * @param height height
     * @param map a double array of size width*height containing occupancy info (optional)
     * @param scale the size of each cell (default is 1.0)
     */
    TravSearch(int width, int height, const double *map = 0, double scale = 1.0);
    
    /**
     * Calculates the cost (usually a height) gradient between two vertices.
     * param v1cost cost of "from" vertex  
     * param v2cost cost of "to" vertex  
     * param ecost cost of edge
     */
    virtual double CalcEdgeCost(double v1cost, double v2cost, double elength);
  };
}


#endif
