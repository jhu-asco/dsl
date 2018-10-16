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

template < class VertexDataT, class EdgeDataT = Empty >
class LpAstar : public Search< VertexDataT, EdgeDataT > {
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
  LpAstar(Graph< VertexDataT, EdgeDataT >& graph,
          const Cost< VertexDataT >& cost);

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
  double plan(std::vector< EdgeT* >& path) override;

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
  void setStart(const VertexT& v) override;

  /**
   * Set goal state
   * this also resets the planner
   * @param v goal vertex
   */
  void addGoal(VertexT& v) override;

  /**
   * Change the cost of edge e
   * @param e edge
   * @param cost new cost
   */
  void changeCost(EdgeT& e, double cost) override;

  /**
   * Change the cost of either all incoming or all outgoing edges
   * connected to verte x
   * @param v vertex
   * @param cost new cost
   * @param in when to modify incoming edges or outgoing edges
   */
  void changeCost(VertexT& v, double cost, bool in = true) override;

  bool inGoalSet(const Vertex< VertexDataT, EdgeDataT >& v) const override {
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
  void updateVertex(VertexT& u);
  void computeShortestPath();
  VertexT* minSucc(double* minRhs, const VertexT& v);
  double minPredRhs(const VertexT& u);
  VertexT* minPred(const VertexT& u);
  double goalSetHeur(const VertexT& s);

  double* calculateExtKey(double* key, VertexT& v);
  double* calculateKey(VertexT& v);
  void insert(VertexT& v);
  void insertExt(VertexT& v, double* key);
  void update(VertexT& v);
  void remove(VertexT& v);
  VertexT* top();
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

  Graph< VertexDataT, EdgeDataT >& graph; ///< graph
  const Cost< VertexDataT >& cost;        ///< cost interface

  std::vector< EdgeT* > changed_edges; ///< newly changed edges

  VertexT* start; ///< start state

  std::map< int, VertexT* > goal_set; ///< all vertices

  fibheap_t open_list; ///< fibonacci heap

  double eps; ///< epsilon for cost comparision

  VertexT* last = 0; ///< last state

  VertexT* goal = nullptr; ///< the first goal found state

  friend class Graph< VertexDataT, EdgeDataT >;
};


// these are needed by fibheap
extern int FIBHEAPKEY_SIZE;
extern fibheapkey_t FIBHEAPKEY_MIN;
extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b);

template < class VertexDataT, class EdgeDataT >
LpAstar< VertexDataT, EdgeDataT >::LpAstar(
    Graph< VertexDataT, EdgeDataT >& graph, const Cost< VertexDataT >& cost)
  : graph(graph), cost(cost), start(0), eps(1e-10) {
  open_list = fibheap_new();
}

template < class VertexDataT, class EdgeDataT >
LpAstar< VertexDataT, EdgeDataT >::~LpAstar() {
  fibheap_delete(open_list);
  changed_edges.clear();
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::reset() {
  typename std::map< int, VertexT* >::iterator vi;
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    vi->second->reset();
  }
  fibheap_clear(open_list);
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::changeCost(EdgeT& edge, double cost) {
  if (nearEqual(cost, edge.cost))
    return;
  edge.cost_change = cost - edge.cost;
  edge.cost = cost;
  changed_edges.push_back(&edge);
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::changeCost(VertexT& vertex,
                                                   double cost,
                                                   bool in) {
  typename std::map< int, EdgeT* >::iterator ei;
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

template < class VertexDataT, class EdgeDataT >
double LpAstar< VertexDataT, EdgeDataT >::plan(std::vector< EdgeT* >& path) {
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
void LpAstar< VertexDataT, EdgeDataT >::setStart(const VertexT& s) {
  start = (VertexT*)&s;
  LpAstar< VertexDataT, EdgeDataT >::reset();
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::addGoal(VertexT& goal) {
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

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::updateVertex(VertexT& u) {
  if (&u != start)
    u.rhs = minPredRhs(u);

  if (u.t == VertexT::Label::kOpen) {
    remove(u);
  }

  if (u.g != u.rhs) {
    insert(u);
  }
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::computeShortestPath() {
  VertexT* u;
  VertexT* s;
  EdgeT* edge;

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
        edge = (EdgeT*)ei.second;
        s = edge->to;
        updateVertex(*s);
      }
    } else {
      u->g = DSL_DBL_MAX;

      for (auto &ei : u->out) {
        edge = (EdgeT*)ei.second;
        s = edge->to;
        updateVertex(*s);
      }
      updateVertex(*u);
    }
  }
}

template < class VertexDataT, class EdgeDataT >
int LpAstar< VertexDataT, EdgeDataT >::plan() {
  VertexT* v;
  int count = 1;

  assert(start);

  if (changed_edges.size()) {
    typename std::vector< EdgeT* >::iterator ei;
    for (ei = changed_edges.begin(); ei != changed_edges.end(); ++ei) {
      EdgeT* edge = *ei;
      v = edge->to;

      updateVertex(*v);
    }
    changed_edges.clear();
  }

  computeShortestPath();

  VertexT* cur = goal;
  do {
    VertexT* prev = minPred(*cur);
    if (!prev) {
      break;
    }
    prev->next = cur;
    cur = prev;
    count++;
  } while (cur != start);

  return count;
}

template < class VertexDataT, class EdgeDataT >
double LpAstar< VertexDataT, EdgeDataT >::minPredRhs(const VertexT& u) {
  double minVal = DSL_DBL_MAX;
  for (auto& ei : u.in) {
    EdgeT* edge = (EdgeT*)ei.second;
    VertexT* s_ = edge->from;

    double val = edge->cost + s_->g;

    if (val < minVal) {
      minVal = val;
    }
  }
  return minVal;
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >*
    LpAstar< VertexDataT, EdgeDataT >::minPred(const VertexT& u) {
  double minVal = DSL_DBL_MAX;
  Vertex< VertexDataT, EdgeDataT >* v = 0;
  for (auto& ei : u.in) {
    EdgeT* edge = (EdgeT*)ei.second;
    VertexT* s_ = edge->from;

    double val = edge->cost + s_->g;

    if (val < minVal) {
      minVal = val;
      v = s_;
    }
  }
  return v;
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >*
    LpAstar< VertexDataT, EdgeDataT >::minSucc(double* minRhs,
                                               const VertexT& s) {
  double minVal = DSL_DBL_MAX;
  VertexT* minSucc = NULL;
  VertexT* s_;
  typename std::map< int, EdgeT* >::const_iterator ei;
  EdgeT* edge;

  for (ei = s.out.begin(); ei != s.out.end(); ++ei) {
    edge = (EdgeT*)ei->second;
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

template < class VertexDataT, class EdgeDataT >
double LpAstar< VertexDataT, EdgeDataT >::goalSetHeur(const VertexT& s) {
  double hmin = DSL_DBL_MAX;
  for (auto &gi : goal_set) {
    VertexT* g = gi.second;
    double h = cost.heur(s.data, g->data);
    if (h < hmin)
      hmin = h;
  }
  return hmin;
}

template < class VertexDataT, class EdgeDataT >
double* LpAstar< VertexDataT, EdgeDataT >::calculateExtKey(double* key,
                                                           VertexT& s) {
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

template < class VertexDataT, class EdgeDataT >
double* LpAstar< VertexDataT, EdgeDataT >::calculateKey(VertexT& s) {
  return calculateExtKey(s.key, s);
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::insert(VertexT& s) {
  insertExt(s, calculateExtKey(s.key, s));
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::insertExt(VertexT& s, double* key) {
  s.open_list_node = fibheap_insert(open_list, (void*)key, &s);
  s.t = VertexT::Label::kOpen;
}

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::update(VertexT& s) {
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

template < class VertexDataT, class EdgeDataT >
void LpAstar< VertexDataT, EdgeDataT >::remove(VertexT& s) {
  if (s.t == VertexT::Label::kClosed)
    return;

  s.t = VertexT::Label::kClosed;
  if (s.open_list_node)
    fibheap_delete_node(open_list, s.open_list_node);
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >* LpAstar< VertexDataT, EdgeDataT >::top() {
  return (VertexT*)fibheap_min(open_list);
}

template < class VertexDataT, class EdgeDataT >
double* LpAstar< VertexDataT, EdgeDataT >::topKey() {
  VertexT* s = top();
  if (!s)
    return NULL;
  return s->key;
}
}
