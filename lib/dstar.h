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

template < class Tv, class Te = Empty >
class Dstar : public Search<Tv, Te> {
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
  Dstar(Graph< Tv, Te >& graph, const Cost< Tv >& cost);

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
  double plan(std::vector< Edge< Tv, Te >* >& path);

  /**
   * Set start state
   * @param v start vertex
   */
  void setStart(const Vertex< Tv, Te >& v);

  /**
   * Set goal state
   * this also resets the planner
   * @param v goal vertex
   */
  void addGoal(Vertex< Tv, Te >& v);

  void setExpandCallback(ExpandCallback expand_callback) {
    this->expand_callback = expand_callback;
  }

  /**
   * Change the cost of edge e
   * @param e edge
   * @param cost new cost
   */
  void changeCost(Edge< Tv, Te >& e, double cost);

  /**
   * Change the cost of either all incoming or all outgoing edges
   * connected to verte x
   * @param v vertex
   * @param cost new cost
   * @param in when to modify incoming edges or outgoing edges
   */
  void changeCost(Vertex< Tv, Te >& v, double cost, bool in = true);

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

  bool inGoalSet(const Vertex< Tv, Te> &v) const {
    return (goal_set.find(v.id) != goal_set.end());
  }

  void updateVertex(Vertex< Tv, Te >& u);
  void computeShortestPath();
  Vertex< Tv, Te >* minSucc(double* minRhs, const Vertex< Tv, Te >& v);
  double* calculateExtKey(double* key, Vertex< Tv, Te >& v);
  double* calculateKey(Vertex< Tv, Te >& v);
  void insert(Vertex< Tv, Te >& v);
  void insertExt(Vertex< Tv, Te >& v, double* key);
  void update(Vertex< Tv, Te >& v);
  void remove(Vertex< Tv, Te >& v);
  Vertex< Tv, Te >* top();
  double* topKey();

private:

  ExpandCallback expand_callback;

  Graph< Tv, Te >& graph; ///< graph
  const Cost< Tv >& cost; ///< cost interface

  std::vector< Edge< Tv, Te >* > changed_edges; ///< newly changed edges

  Vertex< Tv, Te >* start; ///< start state
  std::map< int, Vertex< Tv, Te >* > goal_set; ///< all vertices

  //  Vertex< Tv, Te >* goal;  ///< goal state
  Vertex< Tv, Te >* last;  ///< last state

  double km;          ///< km variable
  fibheap_t open_list; ///< fibonacci heap

  double eps; ///< epsilon for cost comparision

  bool dstar_min; ///< whether to use focussed D* -style min extraction: this was
  /// discovered to reduce vertex expansion (false by default)
  bool goal_bias; ///< whether to employ goal bias heuristic: this can speed-up
  /// the search in easier environments (false by default)

  friend class Graph< Tv, Te >;
};

// these are needed by fibheap
extern int FIBHEAPKEY_SIZE;
extern fibheapkey_t FIBHEAPKEY_MIN;
extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b);

template < class Tv, class Te >
Dstar< Tv, Te >::Dstar(Graph< Tv, Te >& graph, const Cost< Tv >& cost)
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

template < class Tv, class Te >
Dstar< Tv, Te >::~Dstar() {
  fibheap_delete(open_list);
  changed_edges.clear();
}

template < class Tv, class Te >
void Dstar< Tv, Te >::reset() {
  typename std::map< int, Vertex< Tv, Te >* >::iterator vi;
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->reset();
  }
  fibheap_clear(open_list);
  km = 0;
  last = 0;
}

template < class Tv, class Te >
void Dstar< Tv, Te >::setStart(const Vertex< Tv, Te >& s) {
  start = (Vertex< Tv, Te >*)&s;
}

template < class Tv, class Te >
void Dstar< Tv, Te >::addGoal(Vertex< Tv, Te >& goal) {
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

template < class Tv, class Te >
void Dstar< Tv, Te >::changeCost(Edge< Tv, Te >& edge, double cost) {
  if (nearEqual(cost, edge.cost))
    return;
  edge.cost_change = cost - edge.cost;
  edge.cost = cost;
  changed_edges.push_back(&edge);
}

template < class Tv, class Te >
void Dstar< Tv, Te >::changeCost(Vertex< Tv, Te >& vertex,
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

template < class Tv, class Te >
void Dstar< Tv, Te >::updateVertex(Vertex< Tv, Te >& u) {
  if (u.g != u.rhs && u.t == Vertex< Tv, Te >::Label::kOpen) {
    update(u);
  } else if (u.g != u.rhs && u.t != Vertex< Tv, Te >::Label::kOpen) {
    insert(u);
  } else if (u.g == u.rhs && u.t == Vertex< Tv, Te >::Label::kOpen) {
    remove(u);
  }
}

template < class Tv, class Te >
void Dstar< Tv, Te >::computeShortestPath() {
  Vertex< Tv, Te >* u;
  Vertex< Tv, Te >* s;
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
        auto edge = (Edge< Tv, Te >*)ei->second;
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
        auto edge = (Edge< Tv, Te >*)ei->second;
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

template < class Tv, class Te >
double Dstar< Tv, Te >::plan(std::vector< Edge< Tv, Te >* >& path) {
  path.clear();
  plan();
  Vertex< Tv, Te >* cur = start;
  double cost = 0;
  do {
    Edge< Tv, Te >* edge = cur->find(*cur->next, false);
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

template < class Tv, class Te >
int Dstar< Tv, Te >::plan() {
  Vertex< Tv, Te >* cur = start;
  Vertex< Tv, Te >* u, *v;
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
    Vertex< Tv, Te >* next = minSucc(0, *cur);
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

template < class Tv, class Te >
Vertex< Tv, Te >* Dstar< Tv, Te >::minSucc(double* minRhs,
                                            const Vertex< Tv, Te >& s) {
  double minVal = DSL_DBL_MAX;
  Vertex< Tv, Te >* minSucc = NULL;
  Vertex< Tv, Te >* s_;
  Edge< Tv, Te >* edge;

  for (auto ei = s.out.begin(); ei != s.out.end(); ++ei) {
    edge = (Edge< Tv, Te >*)ei->second;
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

template < class Tv, class Te >
double* Dstar< Tv, Te >::calculateExtKey(double* key, Vertex< Tv, Te >& s) {
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

template < class Tv, class Te >
double* Dstar< Tv, Te >::calculateKey(Vertex< Tv, Te >& s) {
  return calculateExtKey(s.key, s);
}

template < class Tv, class Te >
void Dstar< Tv, Te >::insert(Vertex< Tv, Te >& s) {
  insertExt(s, calculateExtKey(s.key, s));
}

template < class Tv, class Te >
void Dstar< Tv, Te >::insertExt(Vertex< Tv, Te >& s, double* key) {
  if (dstar_min) {
    s.r = start;
  } else {
    s.open_list_node = fibheap_insert(open_list, (void*)key, &s);
    s.t = Vertex< Tv, Te >::Label::kOpen;
  }
}

template < class Tv, class Te >
void Dstar< Tv, Te >::update(Vertex< Tv, Te >& s) {

  double key[2];
  // assert(s.t == Vertex<Tv, Te>::Label::kOpen);
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

template < class Tv, class Te >
void Dstar< Tv, Te >::remove(Vertex< Tv, Te >& s) {
  if (s.t == Vertex< Tv, Te >::Label::kClosed) {
    return;
  }

  s.t = Vertex< Tv, Te >::Label::kClosed;
  if (s.open_list_node) {
    fibheap_delete_node(open_list, s.open_list_node);
  }
}

template < class Tv, class Te >
Vertex< Tv, Te >* Dstar< Tv, Te >::top() {
  if (dstar_min) {
    Vertex< Tv, Te >* s;
    while ((s = (Vertex< Tv, Te >*)fibheap_min(open_list))) {
      if (s->r != start) {
        update(*s);
      } else {
        return s;
      }
    }
    return NULL;
  } else {
    return (Vertex< Tv, Te >*)fibheap_min(open_list);
  }
}

template < class Tv, class Te >
double* Dstar< Tv, Te >::topKey() {
  Vertex< Tv, Te >* s = top();
  if (!s) {
    return NULL;
  }
  return s->key;
}
}
