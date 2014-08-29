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
#include "graph.h"
#include "cost.h"

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
 * <a href="https://jshare.johnshopkins.edu/mkobila1/projects/dsl/html/hierarchy.html">Class hierarchy</a>
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
  
  class Search
  {
  public:
    
    /**
     * Initialize dsl with a graph and a cost interface
     * @param graph graph (the method will modify the nodes in the graph
     *                     by changing certain search-related parameters
     *                     stored at the nodes, the data or the original
     *                     edge costs would not be changed)
     * @param cost cost interface
     */
    Search(Graph &graph, const Cost &cost);
    
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
    void SetStart(const Vertex &v);
    
    
    /**
     * Set goal state
     * this also resets the planner
     * @param v goal vertex
     */
    void SetGoal(const Vertex &v);
    
    
    /**
     * Change the cost of edge e
     * @param e edge
     * @param cost new cost
     */
    void ChangeCost(Edge &e, double cost);
        
    /**
     * Set epsilon: this is used to compare cell costs
     * @param eps precision (1e-10 by default)
     */
    void SetEps(double eps) { this->eps = eps; }
            
  protected:
    
    void UpdateVertex(Vertex& u);
    void ComputeShortestPath();
    Vertex* MinSucc(double *minRhs, const Vertex& v);
    double* CalculateExtKey(double *key, Vertex& v);
    double* CalculateKey(Vertex& v);
    void Insert(Vertex& v);
    void InsertExt(Vertex& v, double *key);
    void Update(Vertex& v);
    void Remove(Vertex& v);
    Vertex* Pop();
    Vertex* Top();
    double* TopKey();

    /**
     * Compare two doubles for equality
     * @param a first number
     * @param b second number
     * @return \f$|a-b| < eps\f$ 
     */
    bool Eq(double a, double b) const { return fabs(a - b) < eps; }
    

    Graph &graph;                     ///< graph
    const Cost &cost;                 ///< cost interface
   
    std::vector<Edge*> changedEdges;  ///< newly changed edges
    
    Vertex *start;                    ///< start state
    Vertex *goal;                     ///< goal state
    Vertex *last;                     ///< last state

    double km;                        ///< km variable
    fibheap_t openList;               ///< fibonacci heap
    
    double eps;                       ///< epsilon for cost comparision

  private:
    friend class Graph;
  };
}
#endif
  
