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

namespace dsl {
  
  template<typename T>
    class Vertex;

  template<class T>
    class Graph;

  template<class T>
    class Search;
  
  /**
   * A generic edge b/n two vertices. Contains an edge cost 
   * as well as a costChange variable that could be used
   * for path planning/replanning algorithms such as D* (see Dsl)
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  template<typename T>
  class Edge
  {
 
  public:
    /**
     * Initialize an edge using two vertices and a cost;
     * all parameters are optional since an edge does not
     * physically need to connect two vertices or have associated cost
     * @param from from vertext (optional)
     * @param to to vertex (optional)
     * @param cost cost (optional)
     */
    Edge(Vertex<T> *from = 0, Vertex<T> *to = 0, double cost = 0);
    
    virtual ~Edge();

    int id;            ///< edge id (set internally)
    
    Vertex<T> *from;      ///< from vertex
    Vertex<T> *to;        ///< to vertex

    double cost;       ///< cost
    double costChange; ///< change in cost (used internally)

  private:
    static int s_id;   ///< id counter

    friend class Graph<T>;
    friend class Search<T>;

    //    friend std::ostream& operator<<(std::ostream &os, const Edge<T> &e);
    //    friend std::istream& operator>>(std::istream &is, Edge &e);

  };
  
   template<class T>
    int Edge<T>::s_id = 0;
  
   template<class T>  
     Edge<T>::Edge(Vertex<T> *from, Vertex<T> *to, double cost) :
   id(s_id), from(from), to(to), cost(cost), costChange(-INF) {
     ++s_id;
   }    
  
   template<class T>        
     Edge<T>::~Edge() {
   }
   
 
  /**
   * Output the edge to a stream
   * @param os output stream
   * @param e edge
   * @return the output stream
   template<class T>        
     std::ostream& operator<<(std::ostream &os, const Edge<T> &e) {
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
     std::istream& operator>>(std::istream &is, Edge<T> &e) {   
     return is;
   }
    */

}


#endif
