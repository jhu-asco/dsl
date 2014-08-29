// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GRAPH_H
#define GRAPH_H

#include "vertex.h"
#include "edge.h"

namespace dsl {

  class Search;
  
  /**
   * Generic graph data structure using hashtables for
   * storing vertices and edges for quick add/remove ops
   * 
   * Author: Marin Kobilarov -- Copyright (C) 2004 
   */
  class Graph
  {
  public:
    
    /**
     *   Initialize an empty graph
     */
    Graph();
    
    virtual ~Graph();
    
    /**
     * Add vertex. Vertex is added to map of vertices.
     * No other graph elements are modified.
     * @param v vertex
     */
    void AddVertex(Vertex &v);
    
    
    /**
     * Remove vertex v. Incoming and outgoing edges
     * are removed by default but can be kept intact
     * if paramter re is set ot false
     * @param v vertex
     * @param re remove edges? (true by default)
     * @param del delete/free all removed data? (false by default)
     */
    void RemoveVertex(Vertex &v, bool re = true, bool del = false);
    
    
    /**
     * Add edge. If the from/to vertices
     * of the edge are set then they are modified to
     * incude this edge in their outgoing/incoming respectively
     * list of edges
     * @param e edge
     */
    void AddEdge(Edge &e);
    
    
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
     * @param del delete/free all removed data? (false by default, the idea is that
     *            we care more about speed rather than memory; and these will be
     *            deleted when the graph object is deleted anyways)
     */
    void RemoveEdge(Edge &e, bool update = true, bool del = false);

    /**
     * Checks if a vertex exists
     * @param v the vertex
     * @return true if present in the graph
     */
    bool Exists(const Vertex &v) const;
    
    std::map<int, Edge*> edges;         ///< all edges
    std::map<int, Vertex*> vertices;    ///< all vertices

  protected:
    Search *search;                     ///< search operating on this graph

  private:
    friend class Search;
  };
  
}
#endif
  
  
