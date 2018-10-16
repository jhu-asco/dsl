// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include "vertex.h"
#include <limits>

namespace dsl {

#define DSL_DBL_MAX std::numeric_limits< double >::max()
#define DSL_DBL_MIN std::numeric_limits< double >::min()

template < class VertexDataT, class EdgeDataT >
class Vertex;

/// an empty object used as a default user-supplied edge data
struct Empty {};

/**
 * A generic edge b/n two vertices. Contains an edge cost
 * as well as a cost_change variable that could be used
 * for path planning/replanning algorithms such as D* (see Search)
 * Each vertex stores data of type VertexDataT and each edge stores data of type
 *EdgeDataT,
 * edge data is optional and defaults to the simplest data type, i.e. a bool.
 *
 * Author: Marin Kobilarov
 */
template < class VertexDataT, class EdgeDataT = Empty >
struct Edge {
  using VertexT = Vertex< VertexDataT, EdgeDataT >;

  /**
   * Initialize an edge using two vertices and a cost;
   * all parameters are optional since an edge does not
   * physically need to connect two vertices or have associated cost
   * @param from from vertext (optional)
   * @param to to vertex (optional)
   * @param cost cost (optional)
   */
  Edge(VertexT* from = 0, VertexT* to = 0, double cost = 0)  :
  id(s_id), from(from), to(to), cost(cost), cost_change(DSL_DBL_MIN) {
    ++s_id;
  }

  /**
   * Initialize an edge using its user data, two vertices and a cost;
   * some parameters are optional since an edge does not
   * physically need to connect two vertices or have associated cost
   * @param data data
   * @param from from vertext (optional)
   * @param to to vertex (optional)
   * @param cost cost (optional)
   */
  Edge(const EdgeDataT& data,
       VertexT* from = 0,
       VertexT* to = 0,
       double cost = 0)
    : id(s_id),
      data(data),
      from(from),
      to(to),
      cost(cost),
      cost_change(DSL_DBL_MIN) {
    ++s_id;
  }

  virtual ~Edge() = default;

  int id; ///< edge id (set internally at init)

  EdgeDataT data; ///< edge data

  VertexT* from; ///< from vertex
  VertexT* to;   ///< to vertex

  double cost;       ///< cost
  double cost_change; ///< change in cost (used internally)

private:
  static int s_id; ///< id counter
};

template < class VertexDataT, class EdgeDataT >
int Edge< VertexDataT, EdgeDataT >::s_id = 0;
}
