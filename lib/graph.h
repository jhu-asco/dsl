// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRAPH_H
#define DSL_GRAPH_H

#include "vertex.h"
#include "edge.h"

namespace dsl {

  template<class Tv, class Te>
    class Search;
  
  /**
   * Generic graph data structure using hashtables for
   * storing vertices and edges for quick add/remove ops.
   * Each vertex stores data of type Tv and each edge stores data of type Te,
   * edge data is optional and defaults to the simplest data type, i.e. a bool.   
   * 
   * Author: Marin Kobilarov
   */                                           
  template<class Tv, class Te = bool>
    class Graph {
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
    void AddVertex(Vertex<Tv, Te> &v);
    
    
    /**
     * Remove vertex v. Incoming and outgoing edges
     * are removed by default but can be kept intact
     * if paramter re is set ot false
     * @param v vertex
     * @param re remove edges? (true by default)
     * @param del delete/free all removed data? (false by default)
     */
    void RemoveVertex(Vertex<Tv, Te> &v, bool re = true, bool del = false);
    
    
    /**
     * Add edge. If the from/to vertices
     * of the edge are set then they are modified to
     * incude this edge in their outgoing/incoming respectively
     * list of edges
     * @param e edge
     */
    void AddEdge(Edge<Tv, Te> &e);
    
    
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
    void RemoveEdge(Edge<Tv, Te> &e, bool update = true, bool del = false);

    /**
     * Checks if a vertex exists
     * @param v the vertex
     * @return true if present in the graph
     */
    bool Exists(const Vertex<Tv, Te> &v) const;
    
    std::map<int, Edge<Tv, Te>*> edges;         ///< all edges
    std::map<int, Vertex<Tv, Te>*> vertices;    ///< all vertices

  protected:
    Search<Tv, Te> *search;                     ///< search operating on this graph, this is necessary since in D*, we can remove/add edges to the graph, or change their costs, and the search will be updated accordingly though this field, so that next time it is run it will replan quickly

  private:
    friend class Search<Tv, Te>;
  };


  template<class Tv, class Te>
    Graph<Tv, Te>::Graph() : search(0) {
  }
  

  template<class Tv, class Te>
  Graph<Tv, Te>::~Graph() {
  }
  
  template<class Tv, class Te>
    void Graph<Tv, Te>::AddVertex(Vertex<Tv, Te> &v) {
    vertices[v.id] = &v; 
  }
  
  template<class Tv, class Te>
    void Graph<Tv, Te>::RemoveVertex(Vertex<Tv, Te> &v, bool re, bool del) {
    //  cout << "RemoveVretex: <" << endl;
    if (re) {
      typename std::map<int, Edge<Tv, Te>*>::iterator i;
      // remove all incoming edges
      for (i = v.in.begin(); i != v.in.end(); ++i)
        RemoveEdge(*i->second, true, del);
      // remove all outgoing edges
      for (i = v.out.begin(); i != v.out.end(); ++i)
        RemoveEdge(*i->second, true, del);
    }
    
    // remove from list of vertices
    vertices.erase(v.id);
    
    //  cout << "RemoveVretex: ." << endl;
    //  cout << v << endl;
    
    if (search && search->last)
      search->Remove(v);
    
    if (del)
      delete &v;
  }

  template<class Tv, class Te>
    void Graph<Tv, Te>::AddEdge(Edge<Tv, Te> &e) {
    // attach to start node
    if (e.from)
      e.from->out[e.id] = &e;
    // attach to end node
    if (e.to)
      e.to->in[e.id] = &e;
    // add to list of graph edges
    edges[e.id] = &e;
    
    // if there's an active search
    if (search && search->last) {
      double cost = e.cost;
      e.cost = INF;
      search->ChangeCost(e, cost);
    }
  }
  

  template<class Tv, class Te>
    void Graph<Tv, Te>::RemoveEdge(Edge<Tv, Te> &e, bool update, bool del) {
    // if there's an active search
    if (update && search && search->last) {
      Vertex<Tv, Te> *u = e.from;
      Vertex<Tv, Te> *v = e.to;
      
      double cost = e.cost;
      e.costChange = INF - cost;
      e.cost = INF;
      
      //    cout << e.costChange << " " << e.cost << " " << e.cost - e.costChange << endl;
      //    cout << "u->rhs=" << u->rhs << " v->g=" << v->g << " cost=" << cost << endl;
      
      // update the start vertex
      
      
      if (u && u->openListNode && v && v->openListNode) {
        if (search->Eq(u->rhs, cost + v->g) && u != search->goal) {
          search->MinSucc(&u->rhs, *u);        
        }
        search->UpdateVertex(*u);
      }
    }
    
    // remove from list of edges
    edges.erase(e.id);
    // remove from start node
    if (e.from)
      e.from->out.erase(e.id);
    // remove from end node
    if (e.to)
      e.to->in.erase(e.id);
    
    if (del)
      delete &e;
  }
  
  template<class Tv, class Te>
  bool Graph<Tv, Te>::Exists(const Vertex<Tv, Te> &v) const {
    return (vertices.find(v.id) != vertices.end());
  }  
  
}

#endif
  
  
