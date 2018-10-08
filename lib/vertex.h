// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_VERTEX_H
#define DSL_VERTEX_H

#include <map>
#include <float.h>
#include "fibheap.h"
#include <ostream>
#include <istream>
#include "edge.h"

namespace dsl {

template < class Tv, class Te >
class Graph;

template < class Tv, class Te >
class Dstar;

template < class Tv, class Te >
class LpAstar;


/**
 *  Generic graph vertex containing a list of incoming and outgoing edges
 *  as well as information used for graph search algorithms.
 *  Each vertex stores data of type Tv and each edge store data of type Te,
 *  edge data is optional and defaults to the simplest data type, i.e. a bool.
 *
 *  Author: Marin Kobilarov
 */
template < class Tv, class Te = bool>
class Vertex {
public:
  /**
   * Initialize the vertex
   */
  Vertex();

  /**
   * Initialize the vertex with an optional pointer to user
   * data
   * @param data data
   */
  Vertex(const Tv& data);

  virtual ~Vertex();

  /**
   *  Reset to initial conditions (only the id is kept)
   */
  void Reset();

  /**
   * Finds an edge that connects this
   * vertex with a given vertex v
   * @param v vertex v
   * @param in whether to look among incoming or outgoing edges
   * @return the edge or 0 if none
   */
  Edge< Tv, Te >* Find(const Vertex< Tv, Te >& v, bool in) const;

  int id; ///< vertex id (set internally)

  Tv data; ///< vertex data

  bool succExpanded; ///< is the vertex expanded
  bool predExpanded; ///< is the vertex expanded

  std::map< int, Edge< Tv, Te >* > in;  ///< map of incoming edges
  std::map< int, Edge< Tv, Te >* > out; ///< map of outgoing edges

  Vertex< Tv, Te >* next; ///< next state in a path (used for tracing paths)
  Vertex< Tv, Te >* prev; ///< previous state in a path (used for tracing paths)

protected:
  double rhs; ///< dsl g heuristic values (used internally)
  double g;   ///< dsl rhs heuristic values (used internally)

  static const int NEW = 0;    ///< open list label NEW
  static const int OPEN = 1;   ///< open list label OPEN
  static const int CLOSED = 2; ///< open list label CLOSED

  int t;                  ///< label (one of NEW,OPEN,CLOSED)
  fibnode_t openListNode; ///< heap node associated to this vertex
  double key[2];          ///< heap key

  Vertex< Tv, Te >* r; ///< pointer to a state (from focussed D*)

private:
  static int s_id; ///< id counter

  friend class Graph< Tv, Te >;
  friend class Dstar< Tv, Te >;
  friend class LpAstar< Tv, Te >;

};

template < class Tv, class Te >
int Vertex< Tv, Te >::s_id = 0;

template < class Tv, class Te >
Vertex< Tv, Te >::Vertex()
  : id(s_id), succExpanded(false), predExpanded(false) {
  Reset();
  ++s_id;
}

template < class Tv, class Te >
Vertex< Tv, Te >::Vertex(const Tv& data)
  : id(s_id), data(data), succExpanded(false), predExpanded(false) {
  Reset();
  ++s_id;
}

template < class Tv, class Te >
Vertex< Tv, Te >::~Vertex() {
  in.clear();
  out.clear();
}

template < class Tv, class Te >
void Vertex< Tv, Te >::Reset() {
  next = 0;
  prev = 0;
  rhs = g = DSL_DBL_MAX;
  t = Vertex< Tv, Te >::NEW;
  openListNode = 0;
  key[0] = key[1] = DSL_DBL_MAX;

  r = 0;
}

template < class Tv, class Te >
Edge< Tv, Te >* Vertex< Tv, Te >::Find(const Vertex< Tv, Te >& v,
                                       bool in) const {
  typename std::map< int, Edge< Tv, Te >* >::const_iterator it;

  if (in)
    for (it = this->in.begin(); it != this->in.end(); ++it) {
      Edge< Tv, Te >* e = it->second;
      if (e->from == &v)
        return e;
    }
  else
    for (it = this->out.begin(); it != this->out.end(); ++it) {
      Edge< Tv, Te >* e = it->second;
      if (e->to == &v)
        return e;
    }
  return 0;
}
}

#endif
