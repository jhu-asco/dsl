// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_EDGE_H
#define DSL_EDGE_H

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include "vertex.h"
#include <limits>

namespace dsl {

#define DSL_DBL_MAX std::numeric_limits< double >::max()
#define DSL_DBL_MIN std::numeric_limits< double >::min()

template < class Tv, class Te >
class Vertex;

/// an empty object used as a default user-supplied edge data
struct empty {};

/**
 * A generic edge b/n two vertices. Contains an edge cost
 * as well as a costChange variable that could be used
 * for path planning/replanning algorithms such as D* (see Search)
 * Each vertex stores data of type Tv and each edge stores data of type Te,
 * edge data is optional and defaults to the simplest data type, i.e. a bool.
 *
 * Author: Marin Kobilarov
 */
template < class Tv, class Te = empty >
class Edge {
public:
  /**
   * Initialize an edge using two vertices and a cost;
   * all parameters are optional since an edge does not
   * physically need to connect two vertices or have associated cost
   * @param from from vertext (optional)
   * @param to to vertex (optional)
   * @param cost cost (optional)
   */
  Edge(Vertex< Tv, Te >* from = 0, Vertex< Tv, Te >* to = 0, double cost = 0);

  /**
   * Initialize an edge using its user data, two vertices and a cost;
   * some parameters are optional since an edge does not
   * physically need to connect two vertices or have associated cost
   * @param data data
   * @param from from vertext (optional)
   * @param to to vertex (optional)
   * @param cost cost (optional)
   */
  Edge(const Te& data,
       Vertex< Tv, Te >* from = 0,
       Vertex< Tv, Te >* to = 0,
       double cost = 0);

  virtual ~Edge();

  int id; ///< edge id (set internally at init)

  Te data; ///< edge data

  Vertex< Tv, Te >* from; ///< from vertex
  Vertex< Tv, Te >* to;   ///< to vertex

  double cost;       ///< cost
  double costChange; ///< change in cost (used internally)

private:
  static int s_id; ///< id counter
};

template < class Tv, class Te >
int Edge< Tv, Te >::s_id = 0;

template < class Tv, class Te >
Edge< Tv, Te >::Edge(Vertex< Tv, Te >* from, Vertex< Tv, Te >* to, double cost)
  : id(s_id), from(from), to(to), cost(cost), costChange(DSL_DBL_MIN) {
  ++s_id;
}

template < class Tv, class Te >
Edge< Tv, Te >::Edge(const Te& data,
                     Vertex< Tv, Te >* from,
                     Vertex< Tv, Te >* to,
                     double cost)
  : id(s_id),
    data(data),
    from(from),
    to(to),
    cost(cost),
    costChange(DSL_DBL_MIN) {
  ++s_id;
}

template < class Tv, class Te >
Edge< Tv, Te >::~Edge() {}

/**
 * Output the edge to a stream
 * @param os output stream
 * @param e edge
 * @return the output stream
 template<class T>
 std::ostream& operator<<(std::ostream &os, const Edge<Tv, Te> &e) {
 os << e.id << " ";
 if (e.from)
 os << e.from->id << " ";
 else
 os << "-1 ";
 if (e.to)
 os << e.to->id << " ";
 else
 os << "-1 ";
 os << e.cost << " ";
 os << e.costChange;
 return os;
 };
*/

/**
 * Input the edge from a stream
 * @param is input stream
 * @param e edge
 * @return the input stream
 template<class T>
 std::istream& operator>>(std::istream &is, Edge<Tv, Te> &e) {
 return is;
 }
*/
}

#endif
