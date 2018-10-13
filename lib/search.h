// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include "graph.h"
#include <vector>

/*! \mainpage D*Lite
 * \section Documentation
 * \subsection Intro
 *
 * General implementation of the D*-Lite planner. The algorithm
 * finds the shortest path in a directed graph using A* search
 * and has the ability to quickly replan if any edge costs along the
 * path have changed (using dynamic A*, or D*). A typical
 * application is a robotic vehicle with limited sensing radius that needs
 * to optimally reach a goal in a partially known environment.
 * The robot starts to
 * travel towards the goal and as its sensors refine the terrain map
 * the remaining path is efficiently adjusted or replanned using D*.
 * The package provides a general implementation based
 * on an underlying directed graph (graph ops are done using a fibonacci heap
 * for faster key modifications)
 * as well as a grid/lattice-based implementation (derived from the graph-based
 *one) for
 * search in an environment composed of cells of different "traversibility
 *cost."
 *
 * The library is easy to use and extend. Included is a test executable
 * that demonstrates a typical path planning scenario.
 *
 * The original code was written in 2004 and updated in 2015 to support
 *arbitrary n-dimensional grids
 * and generic data types stored vertices and edges, so that vertices and edges
 *can be regarded as "containers".
 *
 * \subsection Build requirements
 *  g++; cmake
 *
 * \subsection Installation
 *
 * - Get the source and cd into main directory
 * - To compile:
 * - $mkdir build; cd build; cmake ..; make
 * - To test (from main projec directory)
 * -- 1) 2d point test
 * -- $build/bin/test2d bin/map.ppm (look at the generated ppm images to view the
 *result)
 * -- 2) kinematic car test
 * -- $build/bin/cartest bin/lanes.cfg (look at the generated ppm images to view the
 *result)
 *
 * \subsection Class Reference
 * <a href="../../docs/html/hierarchy.html">Class hierarchy</a>
 *
 * \subsection Usage
 *  The underlying structure is a regular directed graph
 *  of vertices and edges that can be added and removed during operation
 *  This implementation serves mostly as a base class for specific
 *  type of problems. Thus it does not define a "real distance" and
 *  "heuristic distance" functions between vertices but allows the user to
 *  supply a cost interface which defines them.
 *
 *  The implementation follows the D*-Lite paper by S.Koenig and M. Likhachev
 *  with several optimizations:
 *
 *  - restructuring of some of the internals allows for reduced
 *    number of heap accesses and edge iterations;
 *  - a fibonacci heap for faster O(1) key decrease
 *  - a "Focussed D*" type of heap extraction is used which
 *      was discovered to be more effective than the current D*-Lite top()
 *      by empirical results this modification results in more than 30% speedup
 *      in time processing (for more complex, i.e. maze-like environments),
 *      results in reduced number of total explored states,
 *      as well as total number of heap accesses. In essense, the gain in
 *      efficiency comes from delaying certain heap operations.
 *      For simple environments there's no siginificant difference
 *  - the ability to dynamically add/remove vertices and edges to enable
 *      anytime and incremental implemnetation
 *
 *
 *   The planner is usually used as follows:
 *   - 0. Create a graph
 *   - 1. Create a cost interface
 *   - 2. Initialize the search using Search(graph, cost)
 *   - 3. Set start vertex using setStart()
 *   - 4. Set goal vertex using addGoal()
 *   - 5. find the optimal path plan()
 *   - 6. follow the generated path until some changes in the graph are observed
 *   - 7. changeCost() -- for every changed cost
 *   - 8. setStart() -- to set the current position
 *   - 9. goto 5 to replan path
 *
 *
 * \subsection Example
 *  see directory test
 * \subsection Author
 *  Copyright (C) 2004, 2015 Marin Kobilarov
 * \subsection Keywords
 * D*, D*-Lite, D-star, "D star", "A*", "D* Lite", "Heuristic Search", "grid
 *planning", "Life-long Planning A*"
 */

namespace dsl {

enum class Method { kDstar, kLpAstar};

template < class Tv, class Te = Empty >
class Search {
public:

  using ExpandCallback = std::function<bool(Vertex< Tv, Te >& from, bool fwd)>;

   virtual ~Search() = default;

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * Internally calls plan()
   * The generated path is set in the provided vector
   * @param path the optimal path
   * @return total cost
   */
   virtual double plan(std::vector< Edge< Tv, Te >* >& path) = 0;

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * The generated path can be obtained by following the next
   * pointers from start vertex until goal vertex is reached
   * This function is implemented by the specific method used
   * @return total number of vertices along the path
   */
   virtual int plan() = 0;

  /**
   * Set start state
   * This function is implemented by the specific method used
   * @param v start vertex
   */
  virtual void setStart(const Vertex< Tv, Te >& v) = 0;

  /**
   * Set goal state
   * this also resets the planner
   * This function is implemented by the specific method used
   * @param v goal vertex
   */
  virtual void addGoal(Vertex< Tv, Te >& v) = 0;

  virtual void changeCost(Edge< Tv, Te >& edge, double cost) = 0;

  virtual void changeCost(Vertex< Tv, Te >& vertex, double cost, bool in) = 0;

  // function invoked for every new expanded node during the exploration
  virtual void setExpandCallback(ExpandCallback expand_callback) = 0;

  virtual void remove(Vertex< Tv, Te >& v) = 0;

  virtual bool inGoalSet(const Vertex< Tv, Te >& v) const = 0;

  virtual bool nearEqual(double a, double b) const = 0;

  virtual Vertex< Tv, Te >* minSucc(double* minRhs, const Vertex< Tv, Te >& v) = 0;

  virtual void updateVertex(Vertex< Tv, Te >& u) = 0;

  // optional pointer to the last start state (only applicable to D*)
  // Since D* replans for a new start, it uses this to updates it heuristic
  virtual Vertex< Tv, Te >* last() { return 0; }

};
}
