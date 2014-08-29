// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDSEARCH3D_H
#define DSL_GRIDSEARCH3D_H

#include "search.h"
#include "gridcost3d.h"

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


namespace dsl {

  /**
   *  Path containing a list of grid points
   */
  class GridPath3D
  {
  public:
  GridPath3D() : pos(0), count(0), len(0) {};
    int *pos;          ///< (x,y,z) point position array of size 3*count
    int count;         ///< number of points along the path 
    double len;        ///< length of path (sum of eucl. distances b/n points)
  };
  

  class GridSearch3D : public Search
  {
 public:
    
    /**
     * The planner allocates all states in the beginning
     * a valid map of size length*width*height should be provided
     * scale is used to multiply each cell value and each neighbor transition value
     * before they are used for planning
     * @param length length
     * @param width width
     * @param height height
     * @param map a double array of size length*width*height containing occupancy info (optional)
     * @param scale the size of each cell (default is 1.0)
     */
    GridSearch3D(int length, int width, int height, const double *map = 0, double scale = 1.0);
    
    
    virtual ~GridSearch3D();
    
    /**
     * Change the cost of an individual cell
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @param cost cost
     */
    void SetCost(int x, int y, int z, double cost);
    
    
    /**
     * Get the cost of an individual cell
     * @param x x
     * @param y y
     * @param z z
     * @return cost
     */
    double GetCost(int x, int y, int z) const;
    
    
    /**
     * Change the costs of all cells at once
     * @param map a length*width*height double array containing occupancy data
     */
    void SetMap(const double *map);
    
    
    /**
     * Set start location in the map
     * @param x x-coord
     * @param y y-coord
     * @param z z-coord
     */
    void SetStart(int x, int y, int z);
    
    
    /**
     * Set goal location in the map
     * @param x x-coord
     * @param y y-coord
     * @param z z-coord
     */
    void SetGoal(int x, int y, int z);
    
    
    /**
     * Compute path b/n start and goal vertices
     * these vertices should be already set
     * @param path the resulting path
     */
    void Plan(GridPath3D &path);
    
    /**
     * Optimize the path produced by Plan
     * @param path an initial obstacle-free path (e.g. one produced by calling Plan(path))
     * @param optPath a new path resulting from optimizing the original path path (i.e. a producing a shorter 
     * path with minimal number of obstacle-free segments)
     */
    void OptPath(const GridPath3D &path, GridPath3D &optPath) const;

    /**
     * Useful method to get the graph vertex at position (x,y)
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @return corresponding vertex or 0 if none there
     */
    Vertex* GetVertex(int x, int y, int z) const; 

    /**
     * Useful method to remove a vertex at (x,y)
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     */
    void RemoveVertex(int x, int y, int z);

    /**
     * Useful method for adding edges b/n vertices
     * @param x1 from x-coordiante
     * @param y1 from y-coordiante
     * @param z1 from z-coordiante
     * @param x2 to x-coordiante
     * @param y2 to y-coordiante
     * @param z2 to z-coordiante
     */    
    void AddEdge(int x1, int y1, int z1, int x2, int y2, int z2);
    
    int length;                        ///< length of cost grid
    int width;                         ///< width of cost grid
    int height;                        ///< height of cost grid
    double scale;                      ///< scale of each cell cost and each transition cost
    double *map;                       ///< cost grid array of size width*height
    
  protected:
    Graph graph;
    GridCost3D cost;
  
    Vertex **vertexMap;                ///< vertex grid array of size width*height
  };
}


#endif
