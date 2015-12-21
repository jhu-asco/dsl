// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDSEARCH_H
#define DSL_GRIDSEARCH_H

#include "search.h"
#include "gridcost.h"
#include "gridconnectivity.h"
#include "gridpath.h"
#include "grid.h"
#include "spline.h"
#include <vector>
#include <iostream>

/**
 *  Grid based D*-Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base Search class
 *
 *  Supports n-dimensional grid with abstract "connectivity" interface
 *
 *  For example in 2D, the cell transition costs can be encoded as the euclidean distance
 *  b/n the centers of the cells (i.e. each cell has 8 neighbors:
 *  the neighbors at N,S,W,E have transition cost of 1, and the 
 *  neighbors at NE, SE, NW, SW have transition costs of sqrt(2)
 *  These transition costs are added to the maximum of the values of
 *  two neighboring cells to compute the cost of their connecting edge.
 *
 *
 *  The planner is used as follows:
 *
 *  1) GridSearch(map, width, height)
 *  2) SetStart(Vector2d(x, y))
 *  3) SetGoal(Vector2d(x, y))
 *  4) Plan(path) and OptPath(path, optPath) -- to plan a path and optimize it
 *  5) follow path until map changes are observed
 *  6) for each change: SetCost(x, y, cost) or change the whole map: SetMap(map)
 *  7) SetStart(Vector2d(x,y)) -- change the start to the current robot position
 *  8) goto 4
 * 
 *
 *  Author: Marin Kobilarov
 */


namespace dsl {

  using namespace Eigen;
  using namespace std;
  
  template<int n>
    class GridSearch : public Search<Cell<n>, GridPath<n> > {
  public:

    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<int, n, 1> Vectorni;
    typedef Vertex<Cell<n>, GridPath<n> > CellVertex;
    typedef Edge<Cell<n>, GridPath<n> > CellEdge;
  
    /**
     * The planner requires a grid, its connectivity, and a cost interface
     * @param grid grid
     * @param connnectivity connectivity interface
     * @param cost cost interface
     * @param expand whether to construct/expand the whole graph at init
     */
    GridSearch(const Grid<n>& grid,  
               const GridConnectivity<n>& connectivity,
               const GridCost<n>& cost,
               bool expand = true);
    
    virtual ~GridSearch();
    
    /**
     * Change the cost of an individual cell
     * @param x position
     * @param cost cost
     */
    bool SetCost(const Vectornd &x, double cost);
    
    /**
     * Get the cost of an individual cell
     * @param x position
     * @return cost
     */
    double GetCost(const Vectornd &x) const;
        
    /**
     * Change the costs of all cells at once
     * @param map a width*height double array containing occupancy data
     */
    //    void SetMap(const double *map);
    
    
    /**
     * Set start location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetStart(const Vectornd &x);
    
    /**
     * Set goal location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetGoal(const Vectornd &x);

    /**
     * Expand successors or predecessors of a given vertex
     * @param from the given vertex
     * @param fwd if true then expand forward in time, i.e. successors, otherwise expand predecessors
     * @return true on success
     */
    bool Expand(CellVertex &from, bool fwd = true);
    
    /**
     * Compute path b/n start and goal vertices
     * these vertices should be already set
     * @param path the resulting path
     * @return true on success
     */
    bool Plan(GridPath<n> &path);
    
    /**
     * Useful method to get the graph vertex at position (x,y)
     * @param x x-coordiante
     * @param y y-coordiante
     * @return corresponding vertex or 0 if none there
     */
    // Vertex<Cell2d>* GetVertex(int x, int y) const; 

    /**
     * Useful method to remove a vertex at (x,y)
     * @param x Euclidean point vector
     */
    bool RemoveCell(const Vectornd &x);

    /**
     * Useful method for adding edges b/n vertices
     * @param x1 from x-coordiante
     * @param y1 from y-coordiante
     * @param x2 to x-coordiante
     * @param y2 to y-coordiante
     */    
    // void AddEdge(int x1, int y1, int x2, int y2);

    /**
     * Experimental path "straightening" function
     * @param path original path
     * @param optPath optimized path
     * @param freeCost (anything above freeCost is considered an obstacle through the path cannoth pass)
     * @param traceStep step with which to trace the path during optimization (should be comparable to the cell size, by defaut is -1 which means the internally the cell size is used)
     */
    void OptPath(const GridPath<n> &path, GridPath<n> &optPath, 
                 double freeCost = 1e-3, double traceStep = -1.0) const;
    /*
    void SplinePath(const GridPath<n> &path, std::vector<Vectornd>& splinePath, 
                 //GridPath<n> &splineCells,
                 double traceStep = 0.1,
                 double scaling = 0.1) const;
    */
  protected:

    Graph<Cell<n>, GridPath<n> > graph;        ///< the underlying graph

    const Grid<n> &grid;                       ///< the grid
    const GridConnectivity<n>& connectivity;   ///< the connectivity interface
    const GridCost<n>& cost;                   ///< the cost interface

    //    bool expand;                   ///< whether to expand all vertices and edges at construction

    CellVertex **vertexMap;              ///< vertex grid array of size width*height
  };


  template<int n>
    GridSearch<n>::GridSearch(const Grid<n>& grid, 
                              const GridConnectivity<n>& connectivity,
                              const GridCost<n>& cost,
                              bool expand) : Search<Cell<n>, GridPath<n> >(graph, cost),
    grid(grid), connectivity(connectivity), cost(cost) {

    vertexMap = new CellVertex*[grid.nc];
    memset(vertexMap, 0, grid.nc*sizeof(CellVertex*));

    if (expand) {
      for (int i = 0; i < grid.nc; ++i) {
        if (grid.cells[i]) {               
          vertexMap[i] = new CellVertex(*grid.cells[i]);
          graph.AddVertex(*vertexMap[i]);        
        }
      }
      
      for (int i = 0; i < grid.nc; ++i) {
        if (grid.cells[i]) {
          
          CellVertex *from = vertexMap[i];
          assert(from);

          std::vector<GridPath<n> > paths;
          connectivity(*grid.cells[i], paths);

          for (int j = 0; j < paths.size(); ++j) {
            const GridPath<n>& path = paths[j];

            int id = grid.Id(path.cells.back().c);
            assert(id >= 0 && id < grid.nc);            
            CellVertex *to = vertexMap[id];
            if (!to) 
              continue;
            
            CellEdge* edge = new CellEdge(path, from, to, path.len);          
            graph.AddEdge(*edge);
          }        
          from->predExpanded = true;
          from->succExpanded = true;
        }
      }
    }
  }


  template<int n>
    bool GridSearch<n>::Expand(CellVertex &from, bool fwd) {
   
    // if this is true then all vertices have already been expanded at construction
    if (fwd && from.succExpanded)
      return true;        

    if (!fwd && from.predExpanded)
      return true;        
   
    // 
    int id = grid.Id(from.data.c);    
    assert(id >= 0 && id < grid.nc);

     std::cout << id << std::endl;

    // cell must exist
    Cell<n> *cell = grid.cells[id];
    assert(cell);
    
    // vertexMap[id] = &from;

    //    if (grid.cells[i]) {               
    // vertexMap[i] = new Vertex<Cell<n> >(*grid.cells[i]);
    // graph.AddVertex(*vertexMap[i]);
      //    }

    std::vector<GridPath<n> > paths;
    connectivity(*cell, paths, fwd);

    for (int j = 0; j < paths.size(); ++j) {
      const GridPath<n> &path = paths[j];
      const Cell<n>& cell = path.cells.back();

      int id = grid.Id(cell.c);
      assert(id >= 0 && id < grid.nc);
            
      //      if (!grid.cells[id])
      //        continue;
      
      // if this vertex doesn't exist, create it and add to graph
      if (!vertexMap[id]) {
        vertexMap[id] = new CellVertex(cell);
        graph.AddVertex(*vertexMap[id]);
      }
      CellVertex *to = vertexMap[id];

      /*
      // fwd: from->to
      Edge<Cell<n>, GridPath<n> >* oldEdge = from.Find(*to, !fwd);
      if (oldEdge) {
        std::cout << "dupl" << std::endl;
        if (oldEdge->cost > costs[j]) {
          graph.RemoveEdge(*oldEdge);
        } else {          
          continue;
        }
      }
      */
      
      // if fwd and incoming edge from->to exists, then do not create a new one
      if (from.Find(*to, !fwd))
        continue;

      CellEdge* edge = fwd ?
        new CellEdge(path, &from, to, path.len) :
        new CellEdge(path, to, &from, path.len);
      graph.AddEdge(*edge);
    }       

    if (fwd)
      from.succExpanded = true;

    if (!fwd)
      from.predExpanded = true;

    return true;
  }  

  template<int n>
    GridSearch<n>::~GridSearch() {
    typename std::map<int, CellEdge*>::iterator ei;
    typename std::map<int, CellVertex*>::iterator vi;
    
    for (ei = graph.edges.begin(); ei != graph.edges.end(); ++ei) {
      delete ei->second;
    }
    for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
      delete vi->second;
    }

    delete[] vertexMap;
  }    

  template<int n>
    bool GridSearch<n>::SetStart(const Vectornd &x) {
    if (!grid.Valid(x)) {
      std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose() << std::endl;
      return false;
    }
    
    int id = grid.Id(x);
    //    std::cout << "id=" << id << std::endl;
    assert(id >= 0 && id < grid.nc);

    // cell must exist
    Cell<n> *cell = grid.cells[id];
    if (!cell) {
      std::cout << "[W] GridSearch:SetStart: cell does not exist x=" << x.transpose() << std::endl;
      return false;
    }    

    // if it's not added previously add it
    if (!vertexMap[id]) {
      vertexMap[id] = new CellVertex(*cell);
      graph.AddVertex(*vertexMap[id]);
    }

    //    Expand(*vertexMap[id]);

    Search<Cell<n>, GridPath<n> >::SetStart(*vertexMap[id]);

    return true;
  }
  
  template<int n>
    bool GridSearch<n>::SetGoal(const Vectornd &x) {
    if (!grid.Valid(x)) {
      std::cout << "[W] GridSearch:SetGoal: invalid x=" << x.transpose() << std::endl;
      return false;
    }
    //  if (x < 0 || x >= width || y < 0 || y >= height)
    //    return;
    int id = grid.Id(x);
    assert(id >= 0 && id < grid.nc);
    
    Cell<n> *cell = grid.cells[id];
    if (!cell) {
      std::cout << "[W] GridSearch:SetGoal: cell does not exist x=" << x.transpose() << std::endl;
      return false;
    }

    // if it's not added previously add it
    if (!vertexMap[id]) {
      vertexMap[id] = new CellVertex(*cell);
      graph.AddVertex(*vertexMap[id]);
    }

    //    Expand(*vertexMap[id], false);

    Search<Cell<n>, GridPath<n> >::SetGoal(*vertexMap[id]); 

    return true;
  }

  
  template<int n>
    bool GridSearch<n>::RemoveCell(const Vectornd &x) {
    int id = grid.Id(x);
    if (id < 0 || id >= grid.nc)
      return false;
    
    CellVertex *v = vertexMap[id];
    if (v) {
      graph.RemoveVertex(*v);
      return true;
    }
    return false;
  }
  
  template<int n>
    bool GridSearch<n>::Plan(GridPath<n>& path) {
    
    path.cells.clear();
    path.len = 0;
    
    std::vector<Edge<Cell<n>, GridPath<n> >* > edgePath;
    //    vector<Edge<Cell<n> , GridPath<n> >* > edgePath;

    Search<Cell<n>, GridPath<n> >::Plan(edgePath);

    //    double len = Search<Cell<n>, GridPath<n> >::Plan(edgePath);

    typename vector<CellEdge*>::iterator it;
    for (it = edgePath.begin(); it != edgePath.end(); ++it) {
      CellEdge *edge = *it;
      path.cells.insert(path.cells.end(), edge->data.cells.begin(), edge->data.cells.end());
      path.len += edge->data.len;
    }
    
    
    return true;
  }
  
  template<int n>
    void GridSearch<n>::OptPath(const GridPath<n> &path, GridPath<n> &optPath, 
                                double freeCost, 
                                double traceStep) const {

    double len = 0;
    
    optPath.cells.clear();
    optPath.len = 0;
    
    if (path.cells.size() == 2) {
      optPath.cells = path.cells;
      optPath.len = path.len;
      return;
    }
    typename vector<Cell<n> >::const_iterator it0;
    it0 = path.cells.begin();
    typename vector<Cell<n> >::const_iterator it1;
    it1 = it0 + 1;

    Vectornd x0 = it0->c;
    Vectornd x1 = it1->c;
    
    Vectornd dx0 = x1 - x0;
    double dn = dx0.norm();
    dx0 /= dn;
    
    optPath.cells.push_back(path.cells[0]);
      
    if (traceStep <= 0)
      traceStep = grid.cs.norm();

    //    for (unsigned int i = 1; i < path.cells.size() - 1; ++i) {
    for (; it1 != path.cells.end(); ++it1) {
      x1 = it1->c;
      Vectornd x2 = (it1 + 1)->c;
      Vectornd dx1 = x2 - x0;
      dn = dx1.norm();
      assert(dn > 1e-12);
      dx1 /= dn;
      
      if ((dx0-dx1).norm() > 1e-16) {
        dn = (x0-x2).norm();
        for (double d = traceStep; d < dn; d += traceStep) {
          Vectornd x = x0 + dx1*d;
          // assert(x > 0);
          // assert(y > 0);
          int id = grid.Id(x); 
          assert(id >= 0 && id < grid.nc);
          //        int yi = min((int)round(y), height - 1);
          //        int xi = min((int)round(x), width - 1);
          if (!grid.cells[id] || grid.cells[id]->cost > freeCost) {
            //	  Cell2d cell((int)x1, (int)y1);
            optPath.cells.push_back(*it1);
            //	  len += sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            x0 = x1;
            break;
          }
        }
        dx0 = x2 - x0;
        dn = dx0.norm();
        dx0 /= dn;
      }
    }
    
    //  Cell2d cell((int)x2, (int)y2);
    optPath.cells.push_back(path.cells.back());
    
    //  len += sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));
    //  optPath.len = len;
  }

  
  //template<int n>
  //  void GridSearch<n>::SplinePath(const GridPath<n> &path, std::vector<Vectornd> &splinePath, 
  //                              /*GridPath<n> &splineCells, */
  //                              double traceStep, double scaling) const {
  //  std::vector<double> steps(path.cells.size()); 
  //  std::vector<std::vector<double> > pts(n, std::vector<double>(path.cells.size()));
  //  for(int i = 0; i < path.cells.size(); i++)
  //  {
  //    steps[i] = i*scaling;
  //    for(int j = 0; j < n; j++)
  //    {
  //      pts[j][i] = path.cells[i].c(j);
  //    }
  //  }

  //  std::vector<Spline<double, double> > sps;
  //  for(int i = 0; i < n; i++)
  //  {
  //    sps.push_back(Spline<double, double>(steps, pts[i]));
  //  }
  //  int count = (path.cells.size()-1)/traceStep;
  //  splinePath.resize(count);
  //  //splineCells.cells.resize(count);
  //  //splineCells.len = 0;

  //  for(int i = 0; i < count; i++)
  //  {
  //    Vectornd pti;
  //    for(int j = 0; j < n; j++)
  //    {
  //      pti(j) = sps[j][i*traceStep*scaling];
  //    }
  //    splinePath[i] = pti;

  //    //std::cout << i*traceStep << std::endl;
  //    //std::cout << pti.transpose() << std::endl;
  //    /*
  //    Cell<n>* cell = grid.Get(splinePath[i]);
  //    if(!cell)
  //    {
  //      std::cout << "dsl::SplinePath: Cell does not exist" << std::endl;
  //      return;
  //    }
  //    splineCells.cells[i] = *(grid.Get(splinePath[i]));
  //    if(i > 0)
  //    {
  //      splineCells.len += (splineCells.cells[i].c-splineCells.cells[i-1].c).norm();
  //    }
  //    */
  //  }
  //}
    

  template<int n>
    double GridSearch<n>::GetCost(const Vectornd &x) const {
    
    int id = grid.Id(x);
    if (id < 0 || id >= grid.nc)
      return 0;
    
    Cell<n> *cell = grid.cells[id];
    if (!cell) {
      std::cout << "[W] GridSearch::SetCost: no cell at position " << x.transpose() << std::endl;
      return 0;
    }
    return cell->cost;
  }

  template<int n>
    bool GridSearch<n>::SetCost(const Vectornd &x, double cost) {

    int id = grid.Id(x);
    if (id < 0 || id >= grid.nc)
      return false;

    Cell<n> *cell = grid.cells[id];
    if (!cell) {
      std::cout << "[W] GridSearch::SetCost: no cell at position " << x.transpose() << std::endl;
      return false;
    }
      
    // if the cost has not changed simply return
    if (Search<Cell<n>, GridPath<n> >::Eq(cost, cell->cost))
      return false;
    
    cell->cost = cost;

    CellVertex *vertex = vertexMap[id];
    assert(vertex);
    
    vertex->data.cost = cost;
    
    // fix all connected edges    
    for (typename std::map<int, CellEdge* >::iterator ein = vertex->in.begin();
         ein != vertex->in.end(); ein++) {
      this->ChangeCost(*ein->second, this->cost.Real(ein->second->from->data, vertex->data));
    }
    
    for (typename std::map<int, CellEdge*>::iterator eout = vertex->out.begin(); 
         eout != vertex->out.end(); eout++) {
      this->ChangeCost(*eout->second, this->cost.Real(vertex->data, eout->second->to->data));
    }
    return true;
  }
}


#endif
