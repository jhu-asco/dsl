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
#include "search.h"
#include "cost.h"
#include "fibheap.h"

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
 *   - 2. Initialize the search using Dstar(graph, cost)
 *   - 3. Set start vertex using setStart()
 *   - 4. Set goal vertex using setGoal()
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
 *planning"
 */

namespace dsl {

template < class VertexDataT, class EdgeDataT = Empty >
class Dstar : public Search< VertexDataT, EdgeDataT > {
public:
  using VertexT = Vertex< VertexDataT, EdgeDataT >;
  using EdgeT = Edge< VertexDataT, EdgeDataT >;

  using ExpandCallback = std::function< bool(VertexT& from, bool fwd) >;

  /**
   * Initialize dsl with a graph and a cost interface
   * @param graph graph (the method will modify the nodes in the graph
   *                     by changing certain search-related parameters
   *                     stored at the nodes, the data or the original
   *                     edge costs would not be changed)
   * @param cost cost interface
   */
  Dstar(Graph< VertexDataT, EdgeDataT >& graph,
        const Cost< VertexDataT >& cost);

  virtual ~Dstar();
  /**
   *  reset the planner to its initial state
   *  this can be used for example if the goal state was changed
   *  but the graph structure is the same
   *  This is usually called internally at init and with every setGoal()
   */
  void reset();

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * The generated path can be obtained by following the next
   * pointers from start vertex until goal vertex is reached
   * @return total number of vertices along the path
   */
  int plan();

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * Internally calls plan()
   * The generated path is set in the provided vector
   * @param path the optimal path
   * @return total cost
   */
  double plan(std::vector< EdgeT* >& path);

  /**
   * Set start state
   * @param v start vertex
   */
  void setStart(const VertexT& v);

  /**
   * Set goal state
   * this also resets the planner
   * @param v goal vertex
   */
  void addGoal(VertexT& v);

  void setExpandCallback(ExpandCallback expand_callback) {
    this->expand_callback = expand_callback;
  }

  /**
   * Change the cost of edge e
   * @param e edge
   * @param cost new cost
   */
  void changeCost(EdgeT& e, double cost);

  /**
   * Change the cost of either all incoming or all outgoing edges
   * connected to verte x
   * @param v vertex
   * @param cost new cost
   * @param in when to modify incoming edges or outgoing edges
   */
  void changeCost(VertexT& v, double cost, bool in = true);

  /**
   * Set epsilon: this is used to compare cell costs
   * @param eps precision (1e-10 by default)
   */
  void setEps(double eps) {
    this->eps = eps;
  }

  /**
   * Number of vertices
   * @return number of vertices
   */
  int vertexCount() const {
    return graph.vertices.size();
  }

  /**
   * Number of edges
   * @return number of edges
   */
  int edgeCount() const {
    return graph.edges.size();
  }

  void setDstarMin(bool dstar_min) { this->dstar_min = dstar_min; }

  void setGoalBias(bool goal_bias) { this->goal_bias = goal_bias; }


  /**
   * Compare two doubles for equality
   * @param a first number
   * @param b second number
   * @return \f$|a-b| < eps\f$
   */
  bool nearEqual(double a, double b) const {
    return std::abs(a - b) < eps;
  }

  bool inGoalSet(const Vertex< VertexDataT, EdgeDataT >& v) const {
    return (goal_set.find(v.id) != goal_set.end());
  }

  void updateVertex(VertexT& u);
  void computeShortestPath();
  VertexT* minSucc(double* minRhs, const VertexT& v);
  double* calculateExtKey(double* key, VertexT& v);
  double* calculateKey(VertexT& v);
  void insert(VertexT& v);
  void insertExt(VertexT& v, double* key);
  void update(VertexT& v);
  void remove(VertexT& v);
  VertexT* top();
  double* topKey();

private:

  ExpandCallback expand_callback;

  Graph< VertexDataT, EdgeDataT >& graph; ///< graph
  const Cost< VertexDataT >& cost;        ///< cost interface

  std::vector< EdgeT* > changed_edges; ///< newly changed edges

  VertexT* start;                     ///< start state
  std::map< int, VertexT* > goal_set; ///< all vertices

  //  VertexT* goal;  ///< goal state
  VertexT* last; ///< last state

  double km;          ///< km variable
  fibheap_t open_list; ///< fibonacci heap

  double eps; ///< epsilon for cost comparision

  bool dstar_min; ///< whether to use focussed D* -style min extraction: this was
  /// discovered to reduce vertex expansion (false by default)
  bool goal_bias; ///< whether to employ goal bias heuristic: this can speed-up
  /// the search in easier environments (false by default)

  friend class Graph< VertexDataT, EdgeDataT >;
};

// these are needed by fibheap
extern int FIBHEAPKEY_SIZE;
extern fibheapkey_t FIBHEAPKEY_MIN;
extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b);

template < class VertexDataT, class EdgeDataT >
Dstar< VertexDataT, EdgeDataT >::Dstar(Graph< VertexDataT, EdgeDataT >& graph,
                                       const Cost< VertexDataT >& cost)
  : graph(graph),
    cost(cost),
    start(0),
    last(0),
    km(0),
    eps(1e-10),
    dstar_min(false),
    goal_bias(false) {
  open_list = fibheap_new();
}

template < class VertexDataT, class EdgeDataT >
Dstar< VertexDataT, EdgeDataT >::~Dstar() {
  fibheap_delete(open_list);
  changed_edges.clear();
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::reset() {
  for (auto vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->reset();
  }
  fibheap_clear(open_list);
  km = 0;
  last = 0;
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::setStart(
    const Vertex< VertexDataT, EdgeDataT >& s) {
  start = (VertexT*)&s;
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::addGoal(VertexT& goal) {
  if (!start) {
    std::cout << "[W] Dstar::setGoal: start should be set first!" << std::endl;
    return;
  }
  // reset planner
  if (!goal_set.size()) {
    reset();
  }

  goal_set[goal.id] = &goal;

  // set goal
  goal.rhs = 0;
  goal.key[0] = cost.heur(start->data, goal.data);
  goal.key[1] = 0;
  insertExt(goal, goal.key);
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::changeCost(EdgeT& edge, double cost) {
  if (nearEqual(cost, edge.cost))
    return;
  edge.cost_change = cost - edge.cost;
  edge.cost = cost;
  changed_edges.push_back(&edge);
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::changeCost(VertexT& vertex,
                                                 double cost,
                                                 bool in) {
  if (in) {
    for (auto ei = vertex.in.begin(); ei != vertex.in.end(); ++ei) {
      changeCost(*ei->second, cost);
    }
  } else {
    for (auto ei = vertex.out.begin(); ei != vertex.out.end(); ++ei) {
      changeCost(*ei->second, cost);
    }
  }
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::updateVertex(VertexT& u) {
  if (u.g != u.rhs && u.t == VertexT::Label::kOpen) {
    update(u);
  } else if (u.g != u.rhs && u.t != VertexT::Label::kOpen) {
    insert(u);
  } else if (u.g == u.rhs && u.t == VertexT::Label::kOpen) {
    remove(u);
  }
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::computeShortestPath() {
  VertexT* u;
  VertexT* s;
  double kold[2];
  double gOld;

  graph.search = this;

  if (fibheap_empty(open_list)) {
    std::cerr << "[W] Dstar::computeShortestPath: open_list is empty -- most "
                 "likely this means that a there is no path between start and "
                 "goal!" << std::endl;
    return;
  }

  while (topKey() && (fibkey_compare(topKey(), calculateKey(*start)) < 0 ||
                      start->rhs != start->g)) {
    u = top();


    if (expand_callback) {
      expand_callback(*u, false);
    }

    kold[0] = u->key[0];
    kold[1] = u->key[1];

    if (fibkey_compare(kold, calculateKey(*u)) < 0) {
      update(*u);
    } else
    if (u->g > u->rhs) {
      u->g = u->rhs;
      remove(*u);

      for (auto ei = u->in.begin(); ei != u->in.end(); ++ei) {
        auto edge = (EdgeT*)ei->second;
        s = edge->from;
        if (!inGoalSet(*s)) {
          s->rhs = std::min(s->rhs, edge->cost + u->g);
        }
        updateVertex(*s);
      }
    } else {
      gOld = u->g;
      u->g = DSL_DBL_MAX;

      for (auto ei = u->in.begin(); ei != u->in.end(); ++ei) {
        auto edge = (EdgeT*)ei->second;
        s = edge->from;

        if (nearEqual(s->rhs, edge->cost + gOld)) {
          if (!inGoalSet(*s)) {
            minSucc(&s->rhs, *s);
          }
        }

        updateVertex(*s);
      }

      if (!inGoalSet(*u))
        minSucc(&u->rhs, *u);
      updateVertex(*u);
    }
  }
}

template < class VertexDataT, class EdgeDataT >
double Dstar< VertexDataT, EdgeDataT >::plan(std::vector< EdgeT* >& path) {
  path.clear();
  plan();
  VertexT* cur = start;
  double cost = 0;
  do {
    EdgeT* edge = cur->find(*cur->next, false);
    if (!edge) {
      path.clear();
      return -1;
    }

    cost += edge->cost;
    path.push_back(edge);
    cur = cur->next;
  } while (!inGoalSet(*cur));

  return cost;
}

template < class VertexDataT, class EdgeDataT >
int Dstar< VertexDataT, EdgeDataT >::plan() {
  VertexT* cur = start;
  VertexT* u, *v;
  int count = 1;

  assert(start);

  if (!last)
    last = start;

  if (changed_edges.size()) {
    km += (cost.heur(last->data, start->data));
    last = start;
    for (auto edge : changed_edges) {
      u = edge->from;
      v = edge->to;

      if (edge->cost_change < 0) {
        if (!inGoalSet(*u))
          // new cost
          u->rhs = std::min(u->rhs, edge->cost + v->g);
      } else {
        // old cost
        if (nearEqual(u->rhs, edge->cost - edge->cost_change + v->g)) {
          if (!inGoalSet(*u)) {
            minSucc(&u->rhs, *u);
          }
        }
      }
      updateVertex(*u);
    }
    changed_edges.clear();
  }

  computeShortestPath();

  do {
    VertexT* next = minSucc(0, *cur);
    cur->next = next;
    if (!next) {
      break;
    }
    next->prev = cur;
    cur = next;
    count++;
  } while (!inGoalSet(*cur));

  return count;
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >*
    Dstar< VertexDataT, EdgeDataT >::minSucc(double* minRhs, const VertexT& s) {
  double minVal = DSL_DBL_MAX;
  VertexT* minSucc = NULL;
  VertexT* s_;
  EdgeT* edge;

  for (auto ei = s.out.begin(); ei != s.out.end(); ++ei) {
    edge = (EdgeT*)ei->second;
    s_ = edge->to;

    double val = edge->cost + s_->g;

    /*
      // TODO(marin) : goal bias can be re-enabled after
      // the functionality below is updated to support multiple goals
       if (goal_bias) {
      if (val < minVal || (val == minVal && minSucc &&
                           cost.real(s_->data, goal->data) <
                               cost.real(minSucc->data, goal->data))) {
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
  if (minRhs) {
    *minRhs = minVal;
  }
  return minSucc;
}

template < class VertexDataT, class EdgeDataT >
double* Dstar< VertexDataT, EdgeDataT >::calculateExtKey(double* key,
                                                         VertexT& s) {
  double m = std::min(s.g, s.rhs);
  if (m == DSL_DBL_MAX) {
    key[0] = DSL_DBL_MAX;
    key[1] = DSL_DBL_MAX;
  } else {
    assert(start);
    key[0] = m + cost.heur(start->data, s.data) + km;
    key[1] = m;
  }

  return key;
}

template < class VertexDataT, class EdgeDataT >
double* Dstar< VertexDataT, EdgeDataT >::calculateKey(VertexT& s) {
  return calculateExtKey(s.key, s);
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::insert(VertexT& s) {
  insertExt(s, calculateExtKey(s.key, s));
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::insertExt(VertexT& s, double* key) {
  if (dstar_min) {
    s.r = start;
  } else {
    s.open_list_node = fibheap_insert(open_list, (void*)key, &s);
    s.t = VertexT::Label::kOpen;
  }
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::update(VertexT& s) {
  double key[2];
  // assert(s.t == Vertex<VertexDataT, EdgeDataT>::Label::kOpen);
  calculateExtKey(key, s);

  if (dstar_min) {
    s.r = start;
  }

  if (fibkey_compare(key, s.key) > 0) {
    fibheap_delete_node(open_list, s.open_list_node);
    s.open_list_node = fibheap_insert(open_list, key, &s);
  } else if (!fibkey_compare(key, s.key)) {
    return;
  } else {
    fibheap_replace_key(open_list, s.open_list_node, key);
  }
  // copy back to s's own key
  // s->open_list_node->key = s->key;
  s.key[0] = key[0];
  s.key[1] = key[1];
}

template < class VertexDataT, class EdgeDataT >
void Dstar< VertexDataT, EdgeDataT >::remove(VertexT& s) {
  if (s.t == VertexT::Label::kClosed) {
    return;
  }

  s.t = VertexT::Label::kClosed;
  if (s.open_list_node) {
    fibheap_delete_node(open_list, s.open_list_node);
  }
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >* Dstar< VertexDataT, EdgeDataT >::top() {
  if (dstar_min) {
    VertexT* s;
    while ((s = (VertexT*)fibheap_min(open_list))) {
      if (s->r != start) {
        update(*s);
      } else {
        return s;
      }
    }
    return NULL;
  } else {
    return (VertexT*)fibheap_min(open_list);
  }
}

template < class VertexDataT, class EdgeDataT >
double* Dstar< VertexDataT, EdgeDataT >::topKey() {
  VertexT* s = top();
  if (!s) {
    return NULL;
  }
  return s->key;
}
}
