// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "vertex.h"
#include "edge.h"

namespace dsl {

template < class Tv, class Te >
class Search;

/**
 * Generic graph data structure using hashtables for
 * storing vertices and edges for quick add/remove ops.
 * Each vertex stores data of type Tv and each edge stores data of type Te,
 * edge data is optional and defaults to empty
 *
 * Author: Marin Kobilarov
 */
template < class Tv, class Te = Empty >
struct Graph {
  /**
   *   Initialize an empty graph
   */
  Graph() = default;

  /**
   * Add vertex. Vertex is added to map of vertices.
   * No other graph elements are modified.
   * @param v vertex
   */
  void addVertex(Vertex< Tv, Te >& v);

  /**
   * Remove vertex v. Incoming and outgoing edges
   * are removed by default but can be kept intact
   * if paramter re is set ot false
   * @param v vertex
   * @param re remove edges? (true by default)
   * @param del delete/free all removed data? (false by default), since
   *        often the data is a pointer to an externally created object
   */
  void removeVertex(Vertex< Tv, Te >& v, bool re = true, bool del = false);

  /**
   * Add edge. If the from/to vertices
   * of the edge are set then they are modified to
   * incude this edge in their outgoing/incoming respectively
   * list of edges
   * @param e edge
   */
  void addEdge(Edge< Tv, Te >& e);

  /**
   * Remove edge.
   * If the from/to vertices
   * of the edge are set then they are modified to
   * remove this edge in their outgoing/incoming respectively
   * list of edges
   * @param e edge
   * @param update whether to update the graph/search data to reflect the
   *        influence of removing this edge if this is done during search.
   *        Normally, one can pass update=false only when cleaning up the
   *        graph upon deletion, for efficiency
   * @param del delete/free all removed data? (false by default, the idea is
   * that
   *            we care more about speed rather than memory; and these will be
   *            deleted when the graph object is deleted anyways)
   */
  void removeEdge(Edge< Tv, Te >& e, bool update = true, bool del = false);

  /**
   * Checks if a vertex exists
   * @param v the vertex
   * @return true if present in the graph
   */
  bool exists(const Vertex< Tv, Te >& v) const;

  std::map< int, Edge< Tv, Te >* > edges;      ///< all edges
  std::map< int, Vertex< Tv, Te >* > vertices; ///< all vertices

  // search operating on this graph, this is necessary
  // since in D*, we can remove/add edges to the graph,
  // or change their costs, and the search will be
  // updated accordingly though this field, so that next
  // time it is run it will replan quickly
  Search< Tv, Te >* search = 0;
};


template < class Tv, class Te >
void Graph< Tv, Te >::addVertex(Vertex< Tv, Te >& v) {
  vertices[v.id] = &v;
}

template < class Tv, class Te >
void Graph< Tv, Te >::removeVertex(Vertex< Tv, Te >& v, bool re, bool del) {
  if (re) {
    typename std::map< int, Edge< Tv, Te >* >::iterator i;
    // remove all incoming edges
    for (i = v.in.begin(); i != v.in.end(); ++i)
      removeEdge(*i->second, true, del);
    // remove all outgoing edges
    for (i = v.out.begin(); i != v.out.end(); ++i)
      removeEdge(*i->second, true, del);
  }

  // remove from list of vertices
  vertices.erase(v.id);

  if (search && search->last())
    search->remove(v);

  if (del)
    delete &v;
}

template < class Tv, class Te >
void Graph< Tv, Te >::addEdge(Edge< Tv, Te >& edge) {
  // attach to start node
  if (edge.from)
    edge.from->out[edge.id] = &edge;
  // attach to end node
  if (edge.to)
    edge.to->in[edge.id] = &edge;
  // add to list of graph edges
  edges[edge.id] = &edge;

  // if there's an active search
  if (search && search->last()) {
    double cost = edge.cost;
    edge.cost = DSL_DBL_MAX;
    search->changeCost(edge, cost);
  }
}

template < class Tv, class Te >
void Graph< Tv, Te >::removeEdge(Edge< Tv, Te >& edge, bool update, bool del) {
  // if there's an active search
  if (update && search && search->last()) {
    Vertex< Tv, Te >* u = edge.from;
    Vertex< Tv, Te >* v = edge.to;

    double cost = edge.cost;
    edge.cost_change = DSL_DBL_MAX - cost;
    edge.cost = DSL_DBL_MAX;

    // update the start vertex
    if (u && u->open_list_node && v && v->open_list_node) {
      if (search->nearEqual(u->rhs, cost + v->g) && !search->inGoalSet(*u)) {
        search->minSucc(&u->rhs, *u);
      }
      search->updateVertex(*u);
    }
  }

  // remove from list of edges
  edges.erase(edge.id);
  // remove from start node
  if (edge.from)
    edge.from->out.erase(edge.id);
  // remove from end node
  if (edge.to)
    edge.to->in.erase(edge.id);

  if (del)
    delete &edge;
}

template < class Tv, class Te >
bool Graph< Tv, Te >::exists(const Vertex< Tv, Te >& vertex) const {
  return (vertices.find(vertex.id) != vertices.end());
}
}
