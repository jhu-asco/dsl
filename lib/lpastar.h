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

/*
Life-long planning A* implementation. This is a forward search with
multiple goals, and supports cost changes.
See search.h for more info.
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
   *  reset the planner to its initial state
   *  this can be used for example if the goal state was changed
   *  but the graph structure is the same
   *  This is usually called internally at init and with every setGoal()
   */
  virtual void reset();

  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * Internally calls plan()
   * The generated path is set in the provided vector
   * @param path the optimal path
   * @return total cost
   */
   double plan(std::vector< Edge< Tv, Te >* >& path) override;


  /**
   * Plan an initial path from start to goal, or if any cost
   * changes are detected since it was last called (any
   * calls to changeCost()) then it replans.
   * The generated path can be obtained by following the next
   * pointers from start vertex until goal vertex is reached
   * @return total number of vertices along the path
   */
  int plan() override;


  void setExpandCallback(ExpandCallback expand_callback) override {
    this->expand_callback = expand_callback;
  }


  /**
   * Set start state
   * @param v start vertex
   */
  void setStart(const Vertex< Tv, Te >& v) override;

  /**
   * Set goal state
   * this also resets the planner
   * @param v goal vertex
   */
  void addGoal(Vertex< Tv, Te >& v) override;

  /**
   * Change the cost of edge e
   * @param e edge
   * @param cost new cost
   */
  void changeCost(Edge< Tv, Te >& e, double cost) override;

  /**
   * Change the cost of either all incoming or all outgoing edges
   * connected to verte x
   * @param v vertex
   * @param cost new cost
   * @param in when to modify incoming edges or outgoing edges
   */
  void changeCost(Vertex< Tv, Te >& v, double cost, bool in = true) override;

  bool inGoalSet(const Vertex< Tv, Te> &v) const override {
    return (goal_set.find(v.id) != goal_set.end());
  }

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

private:

  void updateVertex(Vertex< Tv, Te >& u);
  void computeShortestPath();
  Vertex< Tv, Te >* minSucc(double* minRhs, const Vertex< Tv, Te >& v);
  double minPredRhs(const Vertex< Tv, Te >& u);
  Vertex< Tv, Te >* minPred(const Vertex< Tv, Te >& u);
  double goalSetHeur(const Vertex< Tv, Te >& s);

  double* calculateExtKey(double* key, Vertex< Tv, Te >& v);
  double* calculateKey(Vertex< Tv, Te >& v);
  void insert(Vertex< Tv, Te >& v);
  void insertExt(Vertex< Tv, Te >& v, double* key);
  void update(Vertex< Tv, Te >& v);
  void remove(Vertex< Tv, Te >& v);
  Vertex< Tv, Te >* top();
  double* topKey();

  /**
   * Compare two doubles for equality
   * @param a first number
   * @param b second number
   * @return \f$|a-b| < eps\f$
   */
  bool nearEqual(double a, double b) const override {
    return std::abs(a - b) < eps;
  }

  ExpandCallback expand_callback;

  Graph< Tv, Te >& graph; ///< graph
  const Cost< Tv >& cost; ///< cost interface

  std::vector< Edge< Tv, Te >* > changed_edges; ///< newly changed edges

  Vertex< Tv, Te >* start; ///< start state

  std::map< int, Vertex< Tv, Te >* > goal_set; ///< all vertices

  fibheap_t open_list; ///< fibonacci heap

  double eps; ///< epsilon for cost comparision

  Vertex< Tv, Te >* last = 0;  ///< last state

  Vertex< Tv, Te >* goal = nullptr; ///< the first goal found state

  friend class Graph< Tv, Te >;
};


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
  open_list = fibheap_new();
}

template < class Tv, class Te >
LpAstar< Tv, Te >::~LpAstar() {
  fibheap_delete(open_list);
  changed_edges.clear();
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::reset() {
  typename std::map< int, Vertex< Tv, Te >* >::iterator vi;
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->reset();
  }
  fibheap_clear(open_list);
}


template < class Tv, class Te >
void LpAstar< Tv, Te >::changeCost(Edge< Tv, Te >& edge, double cost) {
  if (nearEqual(cost, edge.cost))
    return;
  edge.cost_change = cost - edge.cost;
  edge.cost = cost;
  changed_edges.push_back(&edge);
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::changeCost(Vertex< Tv, Te >& vertex,
                                  double cost,
                                  bool in) {
  typename std::map< int, Edge< Tv, Te >* >::iterator ei;
  if (in) {
    for (ei = vertex.in.begin(); ei != vertex.in.end(); ++ei) {
      changeCost(*ei->second, cost);
    }
  } else {
    for (ei = vertex.out.begin(); ei != vertex.out.end(); ++ei) {
      changeCost(*ei->second, cost);
    }
  }
}

template < class Tv, class Te >
double LpAstar< Tv, Te >::plan(std::vector< Edge< Tv, Te >* >& path) {
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
void LpAstar< Tv, Te >::setStart(const Vertex< Tv, Te >& s) {
  start = (Vertex< Tv, Te >*)&s;
  LpAstar<Tv, Te>::reset();
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::addGoal(Vertex< Tv, Te >& goal) {

  if (!start) {
    std::cout << "[E] LpAstart::AddGoal: start should be set before goal!" << std::endl;
    return;
  }

  // set goal
  start->rhs = 0;
  if (goal_set.empty()) {
    start->key[0] = cost.heur(start->data, goal.data);
    // insert the once there is at least one goal
    insert(*start);
  } else {
    start->key[0] = std::min(start->key[0], cost.heur(start->data, goal.data));
  }

  goal_set[goal.id] = &goal;
}


template < class Tv, class Te >
void LpAstar< Tv, Te >::updateVertex(Vertex< Tv, Te >& u) {
  if (&u != start)
    u.rhs = minPredRhs(u);

  if (u.t == Vertex< Tv, Te >::Label::kOpen) {
    remove(u);
  }

  if (u.g != u.rhs) {
    insert(u);
  }
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::computeShortestPath() {
  Vertex< Tv, Te >* u;
  Vertex< Tv, Te >* s;
  Edge< Tv, Te >* edge;

  graph.search = this;

  if (fibheap_empty(open_list)) {
    std::cerr << "[W] LpAstar::computeShortestPath: open_list is empty -- most "
                 "likely this means that a there is no path between start and "
                 "goal!" << std::endl;
    return;
  }

  while(1) {

    if (!topKey())
      break;

    bool done = false;
    for (auto &gi : goal_set) {
      goal = gi.second;
      if (fibkey_compare(topKey(), calculateKey(*goal)) >= 0 &&
          goal->rhs == goal->g) {
        done = true;
        break;
      }
    }
    if (done)
      break;

    u = top();
    remove(*u);

    // expand forward
    if (expand_callback) {
      expand_callback(*u, true);
    }

    if (u->g > u->rhs) {
      u->g = u->rhs;

      for (auto& ei : u->out) {
        edge = (Edge< Tv, Te >*)ei.second;
        s = edge->to;
        updateVertex(*s);
      }
    } else {
      u->g = DSL_DBL_MAX;

      for (auto &ei : u->out) {
        edge = (Edge< Tv, Te >*)ei.second;
        s = edge->to;
        updateVertex(*s);
      }
      updateVertex(*u);
    }
  }
}


template < class Tv, class Te >
int LpAstar< Tv, Te >::plan() {
  Vertex< Tv, Te >* v;
  int count = 1;

  assert(start);

  if (changed_edges.size()) {
    typename std::vector< Edge< Tv, Te >* >::iterator ei;
    for (ei = changed_edges.begin(); ei != changed_edges.end(); ++ei) {
      Edge< Tv, Te >* edge = *ei;
      v = edge->to;

      updateVertex(*v);
    }
    changed_edges.clear();
  }

  computeShortestPath();

  Vertex< Tv, Te >* cur = goal;
  do {
    Vertex< Tv, Te >* prev = minPred(*cur);
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
double LpAstar< Tv, Te >::minPredRhs(const Vertex< Tv, Te >& u) {
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
    Vertex<Tv, Te>* LpAstar< Tv, Te >::minPred(const Vertex< Tv, Te >& u) {
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
Vertex< Tv, Te >* LpAstar< Tv, Te >::minSucc(double* minRhs,
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

    if (val < minVal) {
      minVal = val;
      minSucc = s_;
    }
  }
  if (minRhs)
    *minRhs = minVal;
  return minSucc;
}


template < class Tv, class Te >
double LpAstar< Tv, Te >::goalSetHeur(const Vertex< Tv, Te >& s) {
  double hmin = DSL_DBL_MAX;
  for (auto &gi : goal_set) {
    Vertex< Tv, Te >* g = gi.second;
    double h = cost.heur(s.data, g->data);
    if (h < hmin)
      hmin = h;
  }
  return hmin;
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::calculateExtKey(double* key, Vertex< Tv, Te >& s) {
  double m = std::min(s.g, s.rhs);
  if (m == DSL_DBL_MAX) {
    key[0] = DSL_DBL_MAX;
    key[1] = DSL_DBL_MAX;
  } else {
    assert(start);
    key[0] = m + goalSetHeur(s);
    key[1] = m;
  }

  return key;
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::calculateKey(Vertex< Tv, Te >& s) {
  return calculateExtKey(s.key, s);
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::insert(Vertex< Tv, Te >& s) {
  insertExt(s, calculateExtKey(s.key, s));
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::insertExt(Vertex< Tv, Te >& s, double* key) {
  s.open_list_node = fibheap_insert(open_list, (void*)key, &s);
  s.t = Vertex< Tv, Te >::Label::kOpen;
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::update(Vertex< Tv, Te >& s) {

  double key[2];
  calculateExtKey(key, s);

  if (fibkey_compare(key, s.key) > 0) {
    fibheap_delete_node(open_list, s.open_list_node);
    s.open_list_node = fibheap_insert(open_list, key, &s);
  } else if (!fibkey_compare(key, s.key)) {
    return;
  } else {
    fibheap_replace_key(open_list, s.open_list_node, key);
  }
  // copy back to s's own key
  s.key[0] = key[0];
  s.key[1] = key[1];
}

template < class Tv, class Te >
void LpAstar< Tv, Te >::remove(Vertex< Tv, Te >& s) {
  if (s.t == Vertex< Tv, Te >::Label::kClosed)
    return;

  s.t = Vertex< Tv, Te >::Label::kClosed;
  if (s.open_list_node)
    fibheap_delete_node(open_list, s.open_list_node);
}

template < class Tv, class Te >
Vertex< Tv, Te >* LpAstar< Tv, Te >::top() {
  return (Vertex< Tv, Te >*)fibheap_min(open_list);
}

template < class Tv, class Te >
double* LpAstar< Tv, Te >::topKey() {
  Vertex< Tv, Te >* s = top();
  if (!s)
    return NULL;
  return s->key;
}
}
