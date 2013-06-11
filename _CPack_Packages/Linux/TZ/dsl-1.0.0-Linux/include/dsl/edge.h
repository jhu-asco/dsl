#ifndef DSL_EDGE_H
#define DSL_EDGE_H

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>

namespace dsl {
  
  class Vertex;
  
  /**
   * A generic edge b/n two vertices. Contains an edge cost 
   * as well as a costChange variable that could be used
   * for path planning/replanning algorithms such as D* (see Dsl)
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
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
    Edge(Vertex *from = 0, Vertex *to = 0, double cost = 0);
    
    virtual ~Edge();

    int id;            ///< edge id (set internally)
    
    Vertex *from;      ///< from vertex
    Vertex *to;        ///< to vertex

    double cost;       ///< cost
    double costChange; ///< change in cost (used internally)

  private:
    static int s_id;   ///< id counter

    friend class Graph;
    friend class Search;

    friend std::ostream& operator<<(std::ostream &os, const Edge &e);
    friend std::istream& operator>>(std::istream &is, Edge &e);

  };

  /**
   * Output the edge to a stream
   * @param os output stream
   * @param e edge
   * @return the output stream
   */
  std::ostream& operator<<(std::ostream &os, const Edge &e);

  /**
   * Input the edge from a stream
   * @param is input stream
   * @param e edge
   * @return the input stream
   */
  std::istream& operator>>(std::istream &is, Edge &e);
}

#endif
