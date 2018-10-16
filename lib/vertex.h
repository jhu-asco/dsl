// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <float.h>
#include "fibheap.h"
#include <ostream>
#include <istream>
#include "edge.h"

namespace dsl {

template < class VertexDataT, class EdgeDataT >
class Graph;

template < class VertexDataT, class EdgeDataT >
class Dstar;

template < class VertexDataT, class EdgeDataT >
class LpAstar;

/**
 *  Generic graph vertex containing a list of incoming and outgoing edges
 *  as well as information used for graph search algorithms.
 *  Each vertex stores data of type VertexDataT and each edge store data of type
 *EdgeDataT,
 *  edge data is optional and defaults to the simplest data type, i.e. a bool.
 *
 *  Author: Marin Kobilarov
 */
template < class VertexDataT, class EdgeDataT = bool >
struct Vertex {
  enum class Label {kNew, kOpen, kClosed};

  using VertexT = Vertex< VertexDataT, EdgeDataT >;
  using EdgeT = Edge< VertexDataT, EdgeDataT >;

  /**
   * Initialize the vertex
   */
  Vertex();

  /**
   * Initialize the vertex with an optional pointer to user
   * data
   * @param data data
   */
  Vertex(const VertexDataT& data);

  virtual ~Vertex() = default;

  /**
   *  reset to initial conditions (only the id is kept)
   */
  void reset();

  /**
   * finds an edge that connects this
   * vertex with a given vertex v
   * @param v vertex v
   * @param in whether to look among incoming or outgoing edges
   * @return the edge or 0 if none
   */
  EdgeT* find(const VertexT& v, bool in) const;

  int id; ///< vertex id (set internally)

  VertexDataT data; ///< vertex data

  bool succ_expanded; ///< is the vertex expanded
  bool pred_expanded; ///< is the vertex expanded

  std::map< int, EdgeT* > in;  ///< map of incoming edges
  std::map< int, EdgeT* > out; ///< map of outgoing edges

  VertexT* next; ///< next state in a path (used for tracing paths)
  VertexT* prev; ///< previous state in a path (used for tracing paths)

  double rhs; ///< dsl g heuristic values (used internally)
  double g;   ///< dsl rhs heuristic values (used internally)

  Label t;                ///< label
  fibnode_t open_list_node; ///< heap node associated to this vertex
  double key[2];          ///< heap key

  VertexT* r; ///< pointer to a state (from focussed D*)

private:
  static int s_id; ///< id counter

  friend class Graph< VertexDataT, EdgeDataT >;
  friend class Dstar< VertexDataT, EdgeDataT >;
  friend class LpAstar< VertexDataT, EdgeDataT >;
};

template < class VertexDataT, class EdgeDataT >
int Vertex< VertexDataT, EdgeDataT >::s_id = 0;

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >::Vertex()
  : id(s_id), succ_expanded(false), pred_expanded(false) {
  reset();
  ++s_id;
}

template < class VertexDataT, class EdgeDataT >
Vertex< VertexDataT, EdgeDataT >::Vertex(const VertexDataT& data)
  : id(s_id), data(data), succ_expanded(false), pred_expanded(false) {
  reset();
  ++s_id;
}

template < class VertexDataT, class EdgeDataT >
void Vertex< VertexDataT, EdgeDataT >::reset() {
  next = 0;
  prev = 0;
  rhs = g = DSL_DBL_MAX;
  t = Vertex< VertexDataT, EdgeDataT >::Label::kNew;
  open_list_node = 0;
  key[0] = key[1] = DSL_DBL_MAX;
  r = 0;
}

template < class VertexDataT, class EdgeDataT >
Edge< VertexDataT, EdgeDataT >* Vertex< VertexDataT, EdgeDataT >::find(
    const Vertex< VertexDataT, EdgeDataT >& v, bool in) const {
  if (in) {
    for (auto it = this->in.begin(); it != this->in.end(); ++it) {
      EdgeT* e = it->second;
      if (e->from == &v)
        return e;
    }
  } else {
    for (auto it = this->out.begin(); it != this->out.end(); ++it) {
      EdgeT* e = it->second;
      if (e->to == &v)
        return e;
    }
  }
  return 0;
}
}
