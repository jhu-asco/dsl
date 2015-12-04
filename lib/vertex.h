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


  template<class T>
    class Graph;

  template<class T>
    class Search;

  template<class T>
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
  template<class T>
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
    Vertex(const T& data);
    
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
    Edge<T>* Find(const Vertex<T> &v, bool in) const;
    
    int id;                   ///< vertex id (set internally)
    
    T data;                   ///< data
    
    std::map<int, Edge<T>*> in;  ///< map of incoming edges
    std::map<int, Edge<T>*> out; ///< map of outgoing edges
    
    Vertex<T> *next;             ///< next state in a path (used for tracing paths)
    Vertex<T> *prev;             ///< previous state in a path (used for tracing paths)
    
    double rhs;               ///< dsl g heuristic values (used internally)
    double g;                 ///< dsl rhs heuristic values (used internally)

 protected:
    
    static const int NEW = 0;    ///< open list label NEW
    static const int OPEN = 1;   ///< open list label OPEN
    static const int CLOSED = 2; ///< open list label CLOSED    

    int t;                    ///< label (one of NEW,OPEN,CLOSED)
    fibnode_t openListNode;   ///< heap node associated to this vertex
    double key[2];            ///< heap key
    
    Vertex<T> *r;                ///< pointer to a state (from focussed D*)
    
  private:
    static int s_id;          ///< id counter

    friend class Graph<T>;

    friend class Search<T>;
    
    //    friend std::ostream& perator<<(std::ostream &os, const Vertex<T> &v);

    //    friend std::istream& operator>>(std::istream &is, Vertex<T> &v);    
  };

  
  template<class T>
    int Vertex<T>::s_id = 0;
  
  template<class T>
    Vertex<T>::Vertex() : 
  id(s_id) {
    Reset();
    ++s_id;
  }
  
  template<class T>
    Vertex<T>::Vertex(const T& data) : 
  id(s_id), data(data) {
    Reset();
    ++s_id;
  }
  
  template<typename T>
    Vertex<T>::~Vertex() {
    in.clear();
    out.clear();
  }
  
  template<typename T>    
    void Vertex<T>::Reset() {
    next = 0;
    prev = 0;
    rhs = g = INF;
    t = Vertex<T>::NEW;
    openListNode = 0;
    key[0] = key[1] = INF;
    

    r = 0;
  }
  
  template<class T>
    Edge<T>* Vertex<T>::Find(const Vertex<T> &v, bool in) const {

    typename std::map<int, Edge<T>*>::const_iterator it;
    
    if (in)
      for (it = this->in.begin(); it != this->in.end(); ++it) {
        Edge<T> *e = it->second;
        if (e->from == &v)
          return e;
      }
    else
      for (it = this->out.begin(); it != this->out.end(); ++it) {
        Edge<T> *e = it->second;
        if (e->to == &v)
          return e;
      }
    return 0;
  }
  /*
  template<typename T>
  std::ostream& operator<<(std::ostream &os, const Vertex<T> &v)
  {
    os << v.id << " ";
    std::map<int, Edge<T>*>::const_iterator it;
    os << "[";
    for (it = v.in.begin(); it != v.in.end(); ++it) {
      os << it->second->id << " ";
    }
    os << "]<-*->[";
    for (it = v.out.begin(); it != v.out.end(); ++it) {
      os << it->second->id << " ";
    }
    os << "] ";

    if (v.next)
      os << v.next->id << " ";
    else 
      os << "-1 ";
    
    os << "rhs=" << v.rhs << " g=" << v.g << " t=" << v.t << " node=" << v.openListNode << " key=(" << v.key[0] << "," << v.key[1] << ")" << std::endl;;
    

    return os;
  };
  */  

 /**
   * Input the vertex from a stream
   * @param is input stream
   * @param v vertex
   * @return the input stream
  template<typename T>
   std::istream& operator>>(std::istream &is, Vertex<T> &v)
  {   
    return is;
  }
   */


}



#endif
