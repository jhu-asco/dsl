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

///< maximum possible edge cost, treated as infinity
//#define INF DBL_MAX
#define INF (1e16)


namespace dsl {

  class Graph;
  class Edge;
  
  /**
   *  Generic graph vertex containing a list of incoming and outgoing edges
   *  as well as information used for graph search algorithms.
   *  Specific implementations can either extend this class
   *  with extra functionality or simply pass a pointer to 
   *  additional application data.
   *
   *  Marin Kobilarov -- Copyright (C) 2004
   */
  class Vertex
  {
  public:
    
    /**
     * Initialize the vertex with an optional pointer to user
     * data
     * @param data data
     */
    Vertex(void *data = 0);
    
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
    Edge* Find(const Vertex &v, bool in) const;

    int id;                   ///< vertex id (set internally)

    void *data;               ///< application data (optional)
    std::map<int, Edge*> in;  ///< map of incoming edges
    std::map<int, Edge*> out; ///< map of outgoing edges
    
    Vertex *next;             ///< next state in a path (used for tracing paths)
    Vertex *prev;             ///< previous state in a path (used for tracing paths)
    
    double rhs;               ///< dsl g heuristic values (used internally)
    double g;                 ///< dsl rhs heuristic values (used internally)

 protected:
    
    static const int NEW = 0;    ///< open list label NEW
    static const int OPEN = 1;   ///< open list label OPEN
    static const int CLOSED = 2; ///< open list label CLOSED    

    int t;                    ///< label (one of NEW,OPEN,CLOSED)
    fibnode_t openListNode;   ///< heap node associated to this vertex
    double key[2];            ///< heap key
    
#ifdef DSL_DSTAR_MIN
    Vertex *r;                ///< pointer to a state (from focussed D*)
#endif
    
  private:
    static int s_id;          ///< id counter

    friend class Graph;
    friend class Search;
    
    friend std::ostream& operator<<(std::ostream &os, const Vertex &v);
    friend std::istream& operator>>(std::istream &is, Vertex &v);    
  };

  /**
   * Output the vertex to a stream
   * @param os output stream
   * @param v vertex
   * @return the output stream
   */
  std::ostream& operator<<(std::ostream &os, const Vertex &v);

  /**
   * Input the vertex from a stream
   * @param is input stream
   * @param v vertex
   * @return the input stream
   */
  std::istream& operator>>(std::istream &is, Vertex &v);
}

#endif
