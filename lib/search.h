// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_SEARCH_H
#define DSL_SEARCH_H

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include "graph.h"
#include "cost.h"
#include "fibheap.h"

/*! \mainpage D*Lite
 * \section Documentation
 * \subsection Intro
 *
 * General implementation of the D*-Lite planner. The algorithm
 * finds the shortest path in a directed graph using A* search
 * and has the ability to quickly replan if any edge costs along the
 * path have changed (using dynamic A*, or D*). A typical 
 * application is a robotic vehicle with limited sensing radius that needs 
 * to optimally reach a goal in a partially known environment. 
 * The robot starts to
 * travel towards the goal and as its sensors refine the terrain map 
 * the remaining path is efficiently adjusted or replanned using D*.
 * The package provides a general implementation based
 * on an underlying directed graph (graph ops are done using a fibonacci heap
 * for faster key modifications)
 * as well as a grid-based implementation (derived from the graph-based one) for
 * search in an environment composed of cells of different "traversibility cost."
 * 
 * The library is easy to use and extend. Included is a test executable
 * that demonstrates a typical path planning scenario.
 * 
 * \subsection Installation
 * \subsection Build requirements
 *  g++; cmake
 *
 * \subsubsection Download
 *
 * - download: <a href="../../dsl-1.0.0-Source.tar.gz">dsl-1.0.0-Source.tar.gz</a>
 * - To unzip   >: tar xfz dsl-1.0.0-Source.tar.gz
 * - To compile >: cd dsl-1.0.0-Source; mkdir build; cd build; cmake ..; make
 * - To test >: cd test; bin/test ../bin/map.ppm (look at the generated ppm images to view the result)
 *
 * \subsection Class Reference
 * <a href="../../docs/html/hierarchy.html">Class hierarchy</a>
 *
 * \subsection Usage
 *  The underlying structure is a regular directed graph
 *  of vertices and edges that can be added and removed during operation
 *  This implementation serves mostly as a base class for specific
 *  type of problems. Thus it does not define a "real distance" and
 *  "heuristic distance" functions between vertices but allows the user to
 *  supply a cost interface which defines them.
 * 
 *  The implementation follows the D*-Lite paper by S.Koenig and M. Likhachev
 *  with several optimizations: 
 *
 *  - restructuring of some of the internals allows for reduced 
 *    number of heap accesses and edge iterations; 
 *  - a fibonacci heap for faster O(1) key decrease
 *  - a "Focussed D*" type of heap extraction is used which
 *      was discovered to be more effective than the current D*-Lite Top()
 *      by empirical results this modification results in more than 30% speedup
 *      in time processing (for more complex, i.e. maze-like environments), 
 *      results in reduced number of total explored states, 
 *      as well as total number of heap accesses. In essense, the gain in
 *      efficiency comes from delaying certain heap operations.
 *      For simple environments there's no siginificant difference
 *
 *
 *   The planner is usually used as follows:
 *   - 0. Create a graph
 *   - 1. Create a cost interface
 *   - 2. Initialize the search using Search(graph, cost)
 *   - 3. Set start vertex using SetStart()
 *   - 4. Set goal vertex using SetGoal()
 *   - 5. Find the optimal path Plan()
 *   - 6. follow the generated path until some changes in the graph are observed
 *   - 7. ChangeCost() -- for every changed cost
 *   - 8. SetStart() -- to set the current position
 *   - 9. goto 5 to replan path
 *
 *
 * \subsection Example
 *  see directory test
 * \subsection Author
 *  Copyright (C) 2004 Marin Kobilarov 
 * \subsection Keywords
 * D*, D*-Lite, D-star, "D star", "A*", "D* Lite", "Heuristic Search"
 */

namespace dsl {

  template<class T>
    class Search {
  public:
    
    /**
     * Initialize dsl with a graph and a cost interface
     * @param graph graph (the method will modify the nodes in the graph
     *                     by changing certain search-related parameters
     *                     stored at the nodes, the data or the original
     *                     edge costs would not be changed)
     * @param cost cost interface
     */
    Search(Graph<T> &graph, const Cost<T> &cost);
    
    virtual ~Search();
    
    /**
     *  Reset the planner to its initial state
     *  this can be used for example if the goal state was changed
     *  but the graph structure is the same
     *  This is usually called internally at init and with every SetGoal()
     */
    void Reset();
    
    
    /**
     * Plan an initial path from start to goal, or if any cost
     * changes are detected since it was last called (any 
     * calls to ChangeCost()) then it replans.
     * The generated path can be obtained by following the next
     * pointers from start vertex until goal vertex is reached
     * @return total number of vertices along the path
     */
    int Plan();
    
    
    /**
     * Set start state
     * @param v start vertex
     */
    void SetStart(const Vertex<T> &v);
    
    
    /**
     * Set goal state
     * this also resets the planner
     * @param v goal vertex
     */
    void SetGoal(const Vertex<T> &v);
    
    
    /**
     * Change the cost of edge e
     * @param e edge
     * @param cost new cost
     */
    void ChangeCost(Edge<T> &e, double cost);
        
    /**
     * Set epsilon: this is used to compare cell costs
     * @param eps precision (1e-10 by default)
     */
    void SetEps(double eps) { this->eps = eps; }

  protected:
    
    void UpdateVertex(Vertex<T>& u);
    void ComputeShortestPath();
    Vertex<T>* MinSucc(double *minRhs, const Vertex<T>& v);
    double* CalculateExtKey(double *key, Vertex<T>& v);
    double* CalculateKey(Vertex<T>& v);
    void Insert(Vertex<T>& v);
    void InsertExt(Vertex<T>& v, double *key);
    void Update(Vertex<T>& v);
    void Remove(Vertex<T>& v);
    Vertex<T>* Top();
    double* TopKey();

    /**
     * Compare two doubles for equality
     * @param a first number
     * @param b second number
     * @return \f$|a-b| < eps\f$ 
     */
    bool Eq(double a, double b) const { return fabs(a - b) < eps; }
    

    Graph<T> &graph;                     ///< graph
    const Cost<T> &cost;                 ///< cost interface
   
    std::vector<Edge<T>*> changedEdges;  ///< newly changed edges
    
    Vertex<T> *start;                    ///< start state
    Vertex<T> *goal;                     ///< goal state
    Vertex<T> *last;                     ///< last state

    double km;                        ///< km variable
    fibheap_t openList;               ///< fibonacci heap
    
    double eps;                       ///< epsilon for cost comparision

  public:
    bool dstarMin;   ///< whether to use focussed D* -style min extraction: this was discovered to reduce vertex expansion (false by default)
    bool goalBias;   ///< whether to employ goal bias heuristic (false by default)  


  private:
    friend class Graph<T>;
  };


#define DSL_MIN(a,b) ((a<b)?(a):(b))
  
  
  // these are needed by fibheap
  extern int FIBHEAPKEY_SIZE;// = 2*sizeof(double);
  extern fibheapkey_t FIBHEAPKEY_MIN;// = (void*)(double[2]){-INF, -INF};
  
  extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b);
  /*
  {
    assert(a); assert(b);
    double* af = (double *)a;
    double* bf = (double *)b;
    if ((af[0] < bf[0]) || (af[0] == bf[0] && af[1] < bf[1]))
      return -1;
    if (af[0] == bf[0] && af[1] == bf[1])
      return 0;
    return 1;
  }
  */
  
  template<class T>
    Search<T>::Search(Graph<T> &graph, const Cost<T> &cost) : 
  graph(graph),
    cost(cost),
    start(0), 
    goal(0), 
    last(0), 
    km(0), 
    eps(1e-10),
    dstarMin(false),
    goalBias(false) {

    openList = fibheap_new();
  }
  


  template<class T>
    Search<T>::~Search() {
    fibheap_delete(openList);
    changedEdges.clear();
  }



  template<class T>
    void Search<T>::Reset() {
    typename std::map<int,Vertex<T>*>::iterator vi;
    for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
      vi->second->Reset();
    }
    fibheap_clear(openList);
    km = 0;
    last = 0;
  }

  
  template<class T>
    void Search<T>::SetStart(const Vertex<T> &s) {
    start = (Vertex<T>*)&s;
  }
  

  template<class T>
    void Search<T>::SetGoal(const Vertex<T> &s) {
    if (!start) {
      std::cout << "[W] Search::SetGoal: start should be set first!" << std::endl;
      return;
    }
    
    // reset planner
    Reset();
    // set goal
    goal = (Vertex<T>*)&s;
    goal->rhs = 0;
    goal->key[0] = cost.Heur(start->data, goal->data);
    goal->key[1] = 0;
    InsertExt(*goal, goal->key);
  }
  

  template<class T>
    void Search<T>::ChangeCost(Edge<T> &edge, double cost) {
    if (Eq(cost, edge.cost))
      return;
    edge.costChange = cost - edge.cost;
    edge.cost = cost;
    changedEdges.push_back(&edge);
  }
  

  //#define STDOUT_DEBUG
  
  template<class T>
    void Search<T>::UpdateVertex(Vertex<T> &u) {
#ifdef STDOUT_DEBUG
    printf("UpdateVertex: begin\n");
#endif
    if (u.g != u.rhs && u.t == Vertex<T>::OPEN) {
#ifdef STDOUT_DEBUG
      printf("UpdateVertex: 1\n");
#endif
      Update(u);
    } else if (u.g != u.rhs && u.t != Vertex<T>::OPEN) {
#ifdef STDOUT_DEBUG
      printf("UpdateVertex: 2\n");
#endif
      Insert(u);
    } else if (u.g == u.rhs && u.t == Vertex<T>::OPEN) {
#ifdef STDOUT_DEBUG
      printf("UpdateVertex: 3\n");
#endif
      Remove(u);
    }
  }


  template<class T>
    void Search<T>::ComputeShortestPath() {
    Vertex<T> *u;
    Vertex<T> *s;
    double kold[2];
    double gOld;
    typename std::map<int, Edge<T>*>::iterator ei;
    Edge<T>* edge;
    
    graph.search = this;
    
#ifdef STDOUT_DEBUG
    printf("ComputeShortestPath: begin\n");
#endif  
    
    //  assert(!fibheap_empty(openList));
    if (fibheap_empty(openList)) {
      std::cerr << "[W] Search::ComputeShortestPath: openList is empty -- most likely this means that a there is no path between start and goal!" << std::endl;
      return;
    }
    
    while(TopKey() && (fibkey_compare(TopKey(), CalculateKey(*start)) < 0 || start->rhs != start->g)) {
      u = Top();
#ifdef STDOUT_DEBUG
      printf("ComputeShortestPath:Top() -> ");// u->Print(stdout);
#endif
      
      kold[0] = u->key[0];
      kold[1] = u->key[1];
      
      if (fibkey_compare(kold, CalculateKey(*u)) < 0) {
#ifdef STDOUT_DEBUG  
        printf("ComputeShortestPath: 1 -> ");
        printf("old:[%.2f %.2f]\nnew:[%.2f %.2f]\n", kold[0],kold[1], u->key[0], u->key[1]);
#endif
        
        Update(*u);
      } else
      if (u->g > u->rhs) {
#ifdef STDOUT_DEBUG
        printf("ComputeShortestPath: 2\n");
#endif
        
        u->g = u->rhs;
        Remove(*u);
        
        for (ei = u->in.begin(); ei != u->in.end(); ++ei) {
          edge = (Edge<T>*)ei->second;
          s = edge->from;
          if (s != goal)
            s->rhs = DSL_MIN(s->rhs, edge->cost + u->g);
          UpdateVertex(*s);
        }
      } else {
#ifdef STDOUT_DEBUG
        printf("ComputeShortestPath: 3\n");
#endif
        
        gOld = u->g;	
        u->g = INF;
        
        for (ei = u->in.begin(); ei != u->in.end(); ++ei) {
          edge = (Edge<T>*)ei->second;
          s = edge->from;
          if (Eq(s->rhs, edge->cost + gOld))
            if (s != goal)
              MinSucc(&s->rhs, *s);
          UpdateVertex(*s);
        }
        
        if (u != goal)
          MinSucc(&u->rhs, *u);
        UpdateVertex(*u);
      }
    }
  }
  
  
  template<class T>
    int Search<T>::Plan() {
    
    Vertex<T>* cur = start;
    Vertex<T> *u, *v;
    int count = 1;
    
    assert(start);
    
    if (!start->out.size()) {
      std::cerr << "[W] Search::Plan: start vertex has no outgoing edges!" << std::endl;
      return 0;
    }
    
    if (!goal->in.size()) {
      std::cerr << "[W] Search::Plan: goal vertex has no incoming edges!" << std::endl;
      return 0;
    }
    
    if (!last)
      last = start;
    
    if (changedEdges.size()) {
      km += (cost.Heur(last->data, start->data));
      last = start;
      typename std::vector<Edge<T>*>::iterator ei;
      for (ei = changedEdges.begin(); ei != changedEdges.end(); ++ei) {
        Edge<T>* edge = *ei;
        u = edge->from;
        v = edge->to;
        
        if (edge->costChange < 0) {
          if (u != goal) 
            // new cost
            u->rhs = DSL_MIN(u->rhs, edge->cost + v->g);
        } else {
          // old cost
          if (Eq(u->rhs, edge->cost - edge->costChange + v->g)) {
            if (u != goal) {
              MinSucc(&u->rhs, *u);
            }
          }
        }
        UpdateVertex(*u); 
      }
      changedEdges.clear();
    }
    
    ComputeShortestPath();
    
    do {
      Vertex<T> *next = MinSucc(0, *cur);
      cur->next = next;
      if (!next) {
        break;
      }
      next->prev = cur;
      cur = next;
      count++;
    } while(cur != goal);
    
    /*
      do {
      cur->next = MinSucc(0, *cur);
    cur = cur->next;
    if (!cur) {
    break;
    }
    count++;
    } while(cur != goal);
    */
    return count;
  }
  
  template<class T>
    Vertex<T>* Search<T>::MinSucc(double *minRhs, const Vertex<T> &s) {
    
    double minVal = INF;
    Vertex<T>* minSucc = NULL;
    Vertex<T>* s_;
    typename std::map<int, Edge<T>*>::const_iterator ei;
    Edge<T> *edge;
    
    for (ei = s.out.begin(); ei != s.out.end(); ++ei) {
      edge = (Edge<T>*)ei->second;
      s_ = edge->to;
      
      double val = edge->cost + s_->g;

      if (goalBias) {
        if  (val < minVal || 
             (val == minVal && minSucc && cost.Real(s_->data, goal->data) < cost.Real(minSucc->data, goal->data))) {
          minVal = val;
          minSucc = s_;        
        } 
      } else {
        if   (val < minVal) {
          minVal = val;
          minSucc = s_;
        }
      }
    }
    if (minRhs)
      *minRhs = minVal;
    return minSucc;
  }
  
  template<class T>
    double* Search<T>::CalculateExtKey(double *key, Vertex<T> &s) {
    
    double m = DSL_MIN(s.g, s.rhs);
    if (m == INF) {
      key[0] = INF;
      key[1] = INF;
    } else {
      assert(start);
      key[0] = m + cost.Heur(start->data, s.data) + km;
      key[1] = m;
    }
    
    return key;
  }
  
  template<class T>
    double* Search<T>::CalculateKey(Vertex<T> &s) {
    return CalculateExtKey(s.key, s);
  }
  
  template<class T>
    void Search<T>::Insert(Vertex<T> &s) {
    InsertExt(s, CalculateExtKey(s.key, s));
  }
  
  template<class T>
    void Search<T>::InsertExt(Vertex<T> &s, double *key) {
    
    if (dstarMin) {
      s.r = start;
    } else {
      s.openListNode = fibheap_insert(openList, (void*)key, &s);
      s.t = Vertex<T>::OPEN;
    }
  }
  
  template<class T>
    void Search<T>::Update(Vertex<T> &s) {
    double key[2];
    assert(s.t == Vertex<T>::OPEN);
    CalculateExtKey(key, s);
    
    if (dstarMin) {
      s.r = start;
    }  
    
    if (fibkey_compare(key, s.key) > 0) {
      fibheap_delete_node(openList, s.openListNode);
      s.openListNode = fibheap_insert(openList, key, &s); 
    } else 
      if (!fibkey_compare(key, s.key)) {
        return;
      } else {
        fibheap_replace_key(openList, s.openListNode, key);
      }
    // copy back to s's own key
    // s->openListNode->key = s->key;
    s.key[0] = key[0]; s.key[1] = key[1];
  }
  
  
  template<class T>
    void Search<T>::Remove(Vertex<T> &s) {
    if (s.t == Vertex<T>::CLOSED)
      return;
    
    s.t = Vertex<T>::CLOSED;
    if (s.openListNode)
      fibheap_delete_node(openList, s.openListNode);
  }
  
  template<class T>
    Vertex<T>* Search<T>::Top()  {
    if (dstarMin) {
      Vertex<T>* s;
      while((s = (Vertex<T>*)fibheap_min(openList))) {
        if (s->r != start)
          Update(*s);
        else
          return s;
      }
      return NULL;
    } else {
      return (Vertex<T>*)fibheap_min(openList);
    }
  }
  
  template<class T>
    double* Search<T>::TopKey() {  
    Vertex<T>* s = Top();
    if (!s)
      return NULL;
    return s->key;
  }
}

#endif

