// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <iostream>
#include "cost.h"
#include "fibheap.h"
#include "search.h"

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
 *      was discovered to be more effective than the current D*-Lite Top()
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
 *   - 3. Set start vertex using SetStart()
 *   - 4. Set goal vertex using SetGoal()
 *   - 5. Find the optimal path Plan()
 *   - 6. follow the generated path until some changes in the graph are observed
 *   - 7. ChangeCost() -- for every changed cost
 *   - 8. SetStart() -- to set the current position
 *   - 9. goto 5 to replan path
 *
 *
 * \subsection Example
 *  see directory test
 * \subsection Author
 *  Copyright (C) 2004, 2015 Marin Kobilarov
 * \subsection Keywords
 * D*, D*-Lite, D-star, "D star", "A*", "D* Lite", "Heuristic LpAstar", "grid
 *planning"
 */

namespace dsl {

template < class Tv, class Te = Empty >
class LpAstar : public Search<Tv, Te> {
public:

  using ExpandCallback = std::function<bool(Vertex< Tv, Te >& from, bool fwd)>;

  /**
   * Initialize dsl with a graph and a cost interface
   * @param graph graph (the method will modify the nodes in the graph
   *                     by changing certain search-related parameters
   *                     stored at the nodes, the data or the original
   *                     edge costs would not be changed)
   * @param cost cost interface
   */
  LpAstar(Graph< Tv, Te >& graph, const Cost< Tv >& cost);

  virtual ~LpAstar();


  /**
   *  Reset the planner to its initial state
   *  this can be used for example if the goal state was changed
   *  but the graph structure is the same
   *  This is usually called internally at init and with every SetGoal()
   */
  virtual void Reset();

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to ChangeCost()) then it replans.
   * Internally calls Plan()
   * The generated path is set in the provided vector
   * @param path the optimal path
   * @return total cost
   */
   double Plan(std::vector< Edge< Tv, Te >* >& path);


  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to ChangeCost()) then it replans.
   * The generated path can be obtained by following the next
   * pointers from start vertex until goal vertex is reached
   * @return total number of vertices along the path
   */
  int Plan();


  virtual void setExpandCallback(ExpandCallback expand_callback) {
    this->expand_callback = expand_callback;
  }


  /**
   * Set start state
   * @param v start vertex
   */
  void SetStart(const Vertex< Tv, Te >& v);

  /**
   * Set goal state
   * this also resets the planner
   * @param v goal vertex
   */
  void AddGoal(Vertex< Tv, Te >& v);

  /**
   * Change the cost of edge e
   * @param e edge
   * @param cost new cost
   */
  void ChangeCost(Edge< Tv, Te >& e, double cost);

  /**
   * Change the cost of either all incoming or all outgoing edges
   * connected to verte x
   * @param v vertex
   * @param cost new cost
   * @param in when to modify incoming edges or outgoing edges
   */
  void ChangeCost(Vertex< Tv, Te >& v, double cost, bool in = true);

  /**
   * Set epsilon: this is used to compare cell costs
   * @param eps precision (1e-10 by default)
   */
  void SetEps(double eps) {
    this->eps = eps;
  }

  /**
   * Number of vertices
   * @return number of vertices
   */
  int Vertices() const {
    return graph.vertices.size();
  }

  /**
   * Number of edges
   * @return number of edges
   */
  int Edges() const {
    return graph.edges.size();
  }

  bool InGoalSet(const Vertex< Tv, Te> &v) const {
    return (goalSet.find(v.id) != goalSet.end());
  }


private:

  void UpdateVertex(Vertex< Tv, Te >& u);
  void ComputeShortestPath();
  Vertex< Tv, Te >* MinSucc(double* minRhs, const Vertex< Tv, Te >& v);
  double MinPredRhs(const Vertex< Tv, Te >& u);
  Vertex< Tv, Te >* MinPred(const Vertex< Tv, Te >& u);
  double GoalSetHeur(const Vertex< Tv, Te >& s);

  double* CalculateExtKey(double* key, Vertex< Tv, Te >& v);
  double* CalculateKey(Vertex< Tv, Te >& v);
  void Insert(Vertex< Tv, Te >& v);
  void InsertExt(Vertex< Tv, Te >& v, double* key);
  void Update(Vertex< Tv, Te >& v);
  void Remove(Vertex< Tv, Te >& v);
  Vertex< Tv, Te >* Top();
  double* TopKey();

  /**
   * Compare two doubles for equality
   * @param a first number
   * @param b second number
   * @return \f$|a-b| < eps\f$
   */
  bool Eq(double a, double b) const {
    return std::abs(a - b) < eps;
  }


  ExpandCallback expand_callback;

  Graph< Tv, Te >& graph; ///< graph
  const Cost< Tv >& cost; ///< cost interface

  std::vector< Edge< Tv, Te >* > changedEdges; ///< newly changed edges

  Vertex< Tv, Te >* start; ///< start state

  std::map< int, Vertex< Tv, Te >* > goalSet; ///< all vertices

  fibheap_t openList; ///< fibonacci heap

  double eps; ///< epsilon for cost comparision

  Vertex< Tv, Te >* last = 0;  ///< last state

  Vertex< Tv, Te >* goal = nullptr; ///< the first goal found state
  friend class Graph< Tv, Te >;
};


#define DSL_MIN(a, b) ((a < b) ? (a) : (b))

// these are needed by fibheap
extern int FIBHEAPKEY_SIZE;
extern fibheapkey_t FIBHEAPKEY_MIN;
extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b);

template < class Tv, class Te >
LpAstar< Tv, Te >::LpAstar(Graph< Tv, Te >& graph, const Cost< Tv >& cost)
  : graph(graph),
    cost(cost),
    start(0),
    eps(1e-10) {
  openList = fibheap_new();
}

template < class Tv, class Te >
LpAstar< Tv, Te >::~LpAstar() {
  fibheap_delete(openList);
  changedEdges.clear();
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::Reset() {
  typename std::map< int, Vertex< Tv, Te >* >::iterator vi;
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->Reset();
  }
  fibheap_clear(openList);
}


template < class Tv, class Te >
void LpAstar< Tv, Te >::ChangeCost(Edge< Tv, Te >& edge, double cost) {
  if (Eq(cost, edge.cost))
    return;
  edge.costChange = cost - edge.cost;
  edge.cost = cost;
  changedEdges.push_back(&edge);
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::ChangeCost(Vertex< Tv, Te >& vertex,
                                  double cost,
                                  bool in) {
  typename std::map< int, Edge< Tv, Te >* >::iterator ei;
  if (in) {
    for (ei = vertex.in.begin(); ei != vertex.in.end(); ++ei) {
      ChangeCost(*ei->second, cost);
    }
  } else {
    for (ei = vertex.out.begin(); ei != vertex.out.end(); ++ei) {
      ChangeCost(*ei->second, cost);
    }
  }
}

template < class Tv, class Te >
double LpAstar< Tv, Te >::Plan(std::vector< Edge< Tv, Te >* >& path) {
  path.clear();
  Plan();
  Vertex< Tv, Te >* cur = start;
  double cost = 0;
  do {
    Edge< Tv, Te >* edge = cur->Find(*cur->next, false);
    if (!edge) {
      path.clear();
      return -1;
    }

    cost += edge->cost;
    path.push_back(edge);
    cur = cur->next;
  } while (!InGoalSet(*cur));

  return cost;
}


template < class Tv, class Te >
void LpAstar< Tv, Te >::SetStart(const Vertex< Tv, Te >& s) {
  start = (Vertex< Tv, Te >*)&s;
  LpAstar<Tv, Te>::Reset();
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::AddGoal(Vertex< Tv, Te >& goal) {
  //  this->goal = &goal;

  if (!start) {
    std::cout << "[E] LpAstart::AddGoal: start should be set before goal!" << std::endl;
    return;
  }

  // set goal
  start->rhs = 0;
  if (goalSet.empty()) {
    start->key[0] = cost.Heur(start->data, goal.data);
    // insert the once there is at least one goal
    Insert(*start);
  } else {
    start->key[0] = std::min(start->key[0], cost.Heur(start->data, goal.data));
  }

  goalSet[goal.id] = &goal;
}


template < class Tv, class Te >
void LpAstar< Tv, Te >::UpdateVertex(Vertex< Tv, Te >& u) {
  if (&u != start)
    u.rhs = MinPredRhs(u);

  if (u.t == Vertex< Tv, Te >::OPEN) {
    Remove(u);
  }

  if (u.g != u.rhs) {
    Insert(u);
  }
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::ComputeShortestPath() {
  Vertex< Tv, Te >* u;
  Vertex< Tv, Te >* s;
  Edge< Tv, Te >* edge;

  //graph.search = this;

  if (fibheap_empty(openList)) {
    std::cerr << "[W] LpAstar::ComputeShortestPath: openList is empty -- most "
                 "likely this means that a there is no path between start and "
                 "goal!" << std::endl;
    return;
  }

  while(1) {

    if (!TopKey())
      break;

    bool done = false;
    for (auto &gi : goalSet) {
      goal = gi.second;
      if (fibkey_compare(TopKey(), CalculateKey(*goal)) >= 0 &&
          goal->rhs == goal->g) {
        done = true;
        break;
      }
    }
    if (done)
      break;

    u = Top();
    Remove(*u);

    // expand forward
    if (expand_callback)
      expand_callback(*u, true);

    if (u->g > u->rhs) {
      u->g = u->rhs;

      for (auto& ei : u->out) {
        edge = (Edge< Tv, Te >*)ei.second;
        s = edge->to;
        UpdateVertex(*s);
      }
    } else {
      u->g = DSL_DBL_MAX;

      for (auto &ei : u->out) {
        edge = (Edge< Tv, Te >*)ei.second;
        s = edge->to;
        UpdateVertex(*s);
      }
      UpdateVertex(*u);
    }
  }
}


template < class Tv, class Te >
int LpAstar< Tv, Te >::Plan() {
  Vertex< Tv, Te >* v;
  int count = 1;

  assert(start);

  if (changedEdges.size()) {
    typename std::vector< Edge< Tv, Te >* >::iterator ei;
    for (ei = changedEdges.begin(); ei != changedEdges.end(); ++ei) {
      Edge< Tv, Te >* edge = *ei;
      v = edge->to;

      UpdateVertex(*v);
    }
    changedEdges.clear();
  }

  ComputeShortestPath();

  Vertex< Tv, Te >* cur = goal;
  do {
    Vertex< Tv, Te >* prev = MinPred(*cur);
    if (!prev) {
      break;
    }
    prev->next = cur;
    cur = prev;
    count++;
  } while (cur != start);

  return count;
}

template < class Tv, class Te >
double LpAstar< Tv, Te >::MinPredRhs(const Vertex< Tv, Te >& u) {
  double minVal = DSL_DBL_MAX;
  for (auto& ei : u.in) {
    Edge< Tv, Te >* edge = (Edge< Tv, Te >*)ei.second;
    Vertex< Tv, Te >* s_ = edge->from;

    double val = edge->cost + s_->g;

    if (val < minVal) {
      minVal = val;
    }
  }
  return minVal;
}


template < class Tv, class Te >
    Vertex<Tv, Te>* LpAstar< Tv, Te >::MinPred(const Vertex< Tv, Te >& u) {
  double minVal = DSL_DBL_MAX;
  Vertex<Tv, Te>* v = 0;
  for (auto& ei : u.in) {
    Edge< Tv, Te >* edge = (Edge< Tv, Te >*)ei.second;
    Vertex< Tv, Te >* s_ = edge->from;

    double val = edge->cost + s_->g;

    if (val < minVal) {
      minVal = val;
      v = s_;
    }
  }
  return v;
}


template < class Tv, class Te >
Vertex< Tv, Te >* LpAstar< Tv, Te >::MinSucc(double* minRhs,
                                            const Vertex< Tv, Te >& s) {
  double minVal = DSL_DBL_MAX;
  Vertex< Tv, Te >* minSucc = NULL;
  Vertex< Tv, Te >* s_;
  typename std::map< int, Edge< Tv, Te >* >::const_iterator ei;
  Edge< Tv, Te >* edge;

  for (ei = s.out.begin(); ei != s.out.end(); ++ei) {
    edge = (Edge< Tv, Te >*)ei->second;
    s_ = edge->to;

    double val = edge->cost + s_->g;

    /*    if (goalBias) {
      if (val < minVal || (val == minVal && minSucc &&
                           cost.Real(s_->data, goal->data) <
                               cost.Real(minSucc->data, goal->data))) {
        minVal = val;
        minSucc = s_;
      }
    } else {
    */
      if (val < minVal) {
        minVal = val;
        minSucc = s_;
      }
      // }
  }
  if (minRhs)
    *minRhs = minVal;
  return minSucc;
}


template < class Tv, class Te >
double LpAstar< Tv, Te >::GoalSetHeur(const Vertex< Tv, Te >& s) {
  double hmin = DSL_DBL_MAX;
  for (auto &gi : goalSet) {
    Vertex< Tv, Te >* g = gi.second;
    double h = cost.Heur(s.data, g->data);
    if (h < hmin)
      hmin = h;
  }
  return hmin;
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::CalculateExtKey(double* key, Vertex< Tv, Te >& s) {
  double m = DSL_MIN(s.g, s.rhs);
  if (m == DSL_DBL_MAX) {
    key[0] = DSL_DBL_MAX;
    key[1] = DSL_DBL_MAX;
  } else {
    assert(start);
    key[0] = m + GoalSetHeur(s);
    key[1] = m;
  }

  return key;
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::CalculateKey(Vertex< Tv, Te >& s) {
  return CalculateExtKey(s.key, s);
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::Insert(Vertex< Tv, Te >& s) {
  InsertExt(s, CalculateExtKey(s.key, s));
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::InsertExt(Vertex< Tv, Te >& s, double* key) {
  s.openListNode = fibheap_insert(openList, (void*)key, &s);
  s.t = Vertex< Tv, Te >::OPEN;
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::Update(Vertex< Tv, Te >& s) {

  double key[2];
  // assert(s.t == Vertex<Tv, Te>::OPEN);
  CalculateExtKey(key, s);

  if (fibkey_compare(key, s.key) > 0) {
    fibheap_delete_node(openList, s.openListNode);
    s.openListNode = fibheap_insert(openList, key, &s);
  } else if (!fibkey_compare(key, s.key)) {
    return;
  } else {
    fibheap_replace_key(openList, s.openListNode, key);
  }
  // copy back to s's own key
  // s->openListNode->key = s->key;
  s.key[0] = key[0];
  s.key[1] = key[1];
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::Remove(Vertex< Tv, Te >& s) {
  if (s.t == Vertex< Tv, Te >::CLOSED)
    return;

  s.t = Vertex< Tv, Te >::CLOSED;
  if (s.openListNode)
    fibheap_delete_node(openList, s.openListNode);
}

template < class Tv, class Te >
Vertex< Tv, Te >* LpAstar< Tv, Te >::Top() {
  return (Vertex< Tv, Te >*)fibheap_min(openList);
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::TopKey() {
  Vertex< Tv, Te >* s = Top();
  if (!s)
    return NULL;
  return s->key;
}
}
