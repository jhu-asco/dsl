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
//#include "lattice.h"
#include <vector>
#include <iostream>

/**
 *  Grid based D*-Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base Search class
 *
 *  Supports n-dimensional grid with abstract "connectivity" interface
 *
 *  For example in 2D, the cell transition costs can be encoded as the euclidean
 *distance
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

using Eigen::Matrix;
using std::vector;
using std::map;
using std::pair;

/**
 * Grid Planner that can compute the optimal path between two given cells on a
 *grid.
 * The path is a sequence of cells. Each cell in the grid has a cost and in
 *addition
 * there is cost of transitioning between two cells.
 *
 */
template < class PointType, class DataType >
  class GridSearch : public Search< Cell<PointType, DataType> , GridPath< PointType, DataType > > {
public:

  using GridVertexData = Cell<PointType, DataType>;
  using GridEdgeData = GridPath<PointType, DataType>;
  
  //  typedef Matrix< double, n, 1 > Vectornd;
  //  typedef Matrix< int, n, 1 > Vectorni;

  using CellVertex = Vertex< GridVertexData, GridEdgeData>;
  using CellEdge = Edge< GridVertexData, GridEdgeData>;

  using TCell = Cell<PointType, DataType>;
  using TGrid = Grid<PointType, DataType>;
    
  //  typedef Vertex< Cell< n, CellData >, GridPath< n, CellData, Tp > > CellVertex;
  //  typedef Edge< Cell< n, CellData >, GridPath< n, CellData, Tp > > CellEdge;

  /**
   * The planner requires a grid, its connectivity, and a cost interface. By
   * default it will expand the
   * underlying graph as the search proceeds. Alternatively, the whole graph can
   * be expanded at initialization
   * which can be very expensive for large enironments, and is ony useful if: 1)
   * start-up time is not an issue,
   * 2) there is enough memory, 3) the environment is very complex (e.g. a
   * maze).
   * @param grid grid
   * @param connectivity connectivity interface
   * @param cost cost interface
   * @param expand whether to construct/expand the whole graph at init
   * @param trackEdgesThroughCell whether to keep a list of edges passing
   * through each cell, this is useful for changing the costs these edges if the
   * cost of that cell changes, or to remove these edges if the cell is removed
   */
  GridSearch(const Grid< PointType, DataType >& grid,
             const GridConnectivity< PointType, DataType >& connectivity,
             const GridCost< PointType, DataType >& cost,
             bool expand = false,
             bool trackEdgesThroughCell = false);

  virtual ~GridSearch();

  /**
   * Change the cost of an individual cell
   * @param x position
   * @param cost cost
   */
  bool SetCost(const PointType& x, double cost);

  /**
   * Get the cost of an individual cell
   * @param x position
   * @return cost
   */
  double GetCost(const PointType& x) const;

  //    void SetMap(const double *map);

  /**
   * Set start location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool SetStart(const PointType& x);

  /**
   * Set goal location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool SetGoal(const PointType& x);

  /**
   * Expand successors or predecessors of a given vertex. This is mostly used
   * internally
   * during search. In some cases can be "pre-called" externally in advance of
   * searching in this area of
   * the grid, for efficiency purposes.
   * @param from the given vertex
   * @param fwd if true then expand forward in time, i.e. successors, otherwise
   * expand predecessors
   * @return true on success
   */
  bool Expand(CellVertex& from, bool fwd = true);

  /**
   * Compute path b/n start and goal vertices
   * these vertices should be already set
   * @param path the resulting path
   * @param removeDuplicateCells remove duplicate consecutive cells (this might
   * appear if two consecutive edges
   * contain cells in a way that the last cell of the first edge overlaps with
   * the first cell of the next edge)
   * @return true on success
   */
  bool Plan(GridPath< PointType, DataType >& path, bool removeDuplicateCells = true);

  /*
  virtual CellVertex* CreateVertex(int id) = 0;
  
  virtual void DeleteVertex(int id) = 0;
  
  virtual CellVertex* GetVertex(int id) = 0;
  */
  
  /**
   * Useful method to get the graph vertex at position x
   * @param x position
   * @return corresponding vertex or 0 if none there
   */
  CellVertex* GetVertex(const PointType& x) const;

  /**
   * Useful method to remove a vertex at position x
   * @param x Euclidean point vector
   */
  bool RemoveCell(const PointType& x);

  /**
   * Useful method for adding edges b/n vertices
   * @param x1 start point
   * @param x2 end point
   */
  bool AddEdge(const PointType& x1, const PointType& x2);

  /**
   * Experimental path "straightening" function
   * @param path original path
   * @param optPath optimized path
   * @param freeCost (anything above freeCost is considered an obstacle through
   * the path cannoth pass)
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is -1 which means the
   * internally the cell size is used)
   */
  void OptPath(const GridPath< PointType, DataType >& path,
               GridPath< PointType, DataType >& optPath,
               double freeCost = 1e-3,
               double traceStep = -1.0) const;
   

  /**
   * Experimental path "smoothing" function
   * @param path original path
   * @param splinePath spline path
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is 0.1)
  void SplinePath(const GridPath< PointType, DataType >& path,
                  std::vector< PointType >& splinePath,
                  // GridPath<n,CellData> &splineCells,
                  double traceStep = 0.1) const;
   */

private:
  Graph< GridVertexData, GridEdgeData> graph; ///< the underlying graph

  const Grid< PointType, DataType >& grid; ///< the grid
  const GridConnectivity< PointType, DataType >&
      connectivity;              ///< the connectivity interface
  const GridCost< PointType, DataType >& cost; ///< the cost interface


  //  virtual CellVertex* GetVertex(int id) const = 0; 
  CellVertex** vertexMap; ///< vertex grid array TODO: add this in lattice

  bool trackEdgesThroughCell; ///< whether to keep a list of edges passing
  /// through each cell, this is useful for changing
  /// the costs these edges if the cost of that cell
  /// changes, or to remove these edges if the cell
  /// is removed
  std::map< int, vector< CellEdge* > >
      edgesThroughCell; ///< map of list of edges for each cell
};


template < class PointType, class DataType >
GridSearch< PointType, DataType>::GridSearch(
    const Grid< PointType, DataType>& grid,
    const GridConnectivity< PointType, DataType >& connectivity,
    const GridCost< PointType, DataType >& cost,
    bool expand,
    bool trackEdgesThroughCell)
    : Search< Cell<PointType, DataType>, GridPath<PointType, DataType > >(graph, cost),
    grid(grid),
    connectivity(connectivity),
    cost(cost),
    trackEdgesThroughCell(trackEdgesThroughCell) {

      vertexMap = new CellVertex* [grid.nc];
      memset(vertexMap, 0, grid.nc * sizeof(CellVertex*));

  if (expand) {
    // initialize and add all vertices
    //    vector<Cell> cells;
    //    grid.GetCells(cells);
    for (int i = 0; i < grid.nc; ++i) {
      const TCell *cell = grid.Get(i);
      if (cell) {
        //      if (connectivity.Free(cell.data)) {
        vertexMap[i] = new CellVertex(*cell);
        graph.AddVertex(*vertexMap[i]);
      }
    }

    // expand the successors of each vertex
    //    for (auto cell : grid.cells) {
    //      if (connectivity.Free(cell.data)) {
    for (int i = 0; i < grid.nc; ++i) {
      const TCell *cell = grid.Get(i);
      if (cell) {
        CellVertex* from = vertexMap[i];
        assert(from);

        // generate successor paths
        std::vector< GridPath< PointType, DataType> > paths;
        connectivity(*cell, paths);

        // add LatticeSearch::GetVertex(PointType &x);
        
        for (auto&& path : paths) {
          // find cell where the end of the path falls
          int id = grid.Id(path.cells.back().c);
          assert(id >= 0 && id < grid.nc);
          CellVertex* to = vertexMap[id];
          if (!to)
            continue;

          /*
          if (path.cells.size() >= 2) {
            path.cost = 0;
            typename vector<Cell<n, CellData> >::iterator cit = path.cells.begin();
            for (; (cit+1) != path.cells.end(); ++cit) {
              path.cost += cost.Real(*cit, *(cit+1));
            }
          } else {
            path.cost = cost.Real(from->data, to->data);
          }
          */

          // path.cost = cost.Real(from->data, to->data);

          CellEdge* edge = new CellEdge(path, from, to, path.cost);
          graph.AddEdge(*edge);

          /*
          if (trackEdgesThroughCell) {
            // iterate through path and add the new edge to the list
            // of edges passing through each cell of the path
            for (auto&& cells : path.cells) {
              int id = grid.Id(cells.first);
              //              typename map< int, vector< CellEdge* > >::iterator etcit;
              auto etcit = edgesThroughCell.find(id);
              if (etcit == edgesThroughCell.end()) {
                vector< CellEdge* > edges;
                edges.push_back(edge);
                edgesThroughCell[id] = edges;
              } else {
                etcit->second.push_back(edge);
              }
            }
          }
          */
        }
        from->predExpanded = true;
        from->succExpanded = true;
      }
    }
  }
}

template < class PointType, class DataType >
    bool GridSearch< PointType, DataType>::Expand(CellVertex& from, bool fwd) {
  // if this is true then all vertices have already been expanded at
  // construction
  if (fwd && from.succExpanded)
    return true;

  if (!fwd && from.predExpanded)
    return true;

  //
  int id = grid.Id(from.data.c);
  assert(id >= 0 && id < grid.nc);

  // cell must exist
  const Cell<PointType, DataType>* cell = grid.Get(id);
  assert(cell);
  
  std::vector< GridPath< PointType, DataType > > paths;
  connectivity(*cell, paths, fwd);

  for (auto&& path : paths) {
    TCell& cell = path.cells.back();

    int id = grid.Id(cell.c);
    assert(id >= 0 && id < grid.nc);

    //      if (!grid.cells[id])
    //        continue;

    // if this vertex doesn't exist, create it and add to graph
    if (!vertexMap[id]) {
      vertexMap[id] = new CellVertex(cell);
      graph.AddVertex(*vertexMap[id]);
    }
    CellVertex* to = vertexMap[id];

    /*
    // fwd: from->to
    Edge<Cell<n,CellData>, GridPath<n,CellData> >* oldEdge = from.Find(*to, !fwd);
    if (oldEdge) {
      std::cout << "dupl" << std::endl;
      if (oldEdge->cost > paths[j].len) {
        graph.RemoveEdge(*oldEdge);
      } else {
        continue;
      }
    }
    */

    // if fwd and incoming edge from->to exists, then do not create a new one
    if (from.Find(*to, !fwd))
      continue;

    /*
    if (path.cells.size() >= 2) {
      path.cost = 0;
      typename vector<Cell<n, CellData> >::iterator cit = path.cells.begin();
      for (; (cit+1) != path.cells.end(); ++cit) {
        path.cost += cost.Real(*cit, *(cit+1));
      }
    } else {
      path.cost = cost.Real(from.data, to->data);
    }
    */

    //  path.cost = cost.Real(from.data, to->data);

    CellEdge* edge = fwd ? new CellEdge(path, &from, to, path.cost) :
                           new CellEdge(path, to, &from, path.cost);
    graph.AddEdge(*edge);

    /*
    if (trackEdgesThroughCell) {
      // iterate through path and add the new edge to the list
      // of edges passing through each cell of the path
      for (auto cell : path.cell) {
        int id = grid.Id(cell.first);
        typename map< int, vector< CellEdge* > >::iterator etcit;
        etcit = edgesThroughCell.find(id);
        if (etcit == edgesThroughCell.end()) {
          vector< CellEdge* > edges;
          edges.push_back(edge);
          edgesThroughCell[id] = edges;
        } else {
          etcit->second.push_back(edge);
        }
      }
    }
    */
  }

  if (fwd)
    from.succExpanded = true;
  else
    from.predExpanded = true;

  return true;
}

template < class PointType, class DataType >
    GridSearch< PointType, DataType >::~GridSearch() {

  for (auto edge : graph.edges) {
    delete edge.second;
  }
  for (auto vertex : graph.vertices) {
    delete vertex.second;
  }

  delete[] vertexMap;
}

template < class PointType, class DataType >
bool GridSearch< PointType, DataType >::SetStart(const PointType& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);
  //    std::cout << "id=" << id << std::endl;
  //  assert(id >= 0 && id < grid.nc);

  // cell must exist
  const Cell<PointType, DataType >* cell = grid.Get(id);
  if (!cell) {
    std::cout << "[W] GridSearch:SetStart: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }
  
  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  //    Expand(*vertexMap[id]);

  Search< Cell<PointType, DataType> , GridPath<PointType, DataType> >::SetStart(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType >
    bool GridSearch< PointType, DataType >::SetGoal(const PointType& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetGoal: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);
  assert(id >= 0 && id < grid.nc);

  const TCell* cell = grid.Get(id);
  assert(cell);
  
  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  Search< Cell<PointType, DataType> , GridPath<PointType, DataType>  >::SetGoal(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType >
    bool GridSearch< PointType, DataType >::RemoveCell(const PointType& x) {
  int id = grid.Id(x);  
  if (id < 0 || id >= grid.nc) {
    std::cout << "[W] GridSearch::RemoveCell: id=" << id << " out of bounds!" << std::endl;    
    return false;
  }
    
  CellVertex* v = vertexMap[id];
  if (v) {
    graph.RemoveVertex(*v);
    delete v;
    vertexMap[id] = 0;
  } else {
    std::cout << "[W] GridSearch::RemoveCell: vertex with id=" << id << " doesnt exist!" << std::endl;
    return false;
  }

  /*
  if (grid.cells[id]) {
    delete grid.cells[id];
    grid.cells[id] = 0;
  } else {
    return false;
  }
  */
  
  /*
  // if edges through cells are tracked, then go through them and remove them
  if (trackEdgesThroughCell) {
    typename map< int, vector< CellEdge* > >::iterator ceit =
        edgesThroughCell.find(id);
    if (ceit != edgesThroughCell.end()) {
      vector< CellEdge* > edges = ceit->second;
      typename vector< CellEdge* >::iterator ei;
      for (ei = edges.begin(); ei != edges.end(); ++ei) {
        this->graph.RemoveEdge(**ei);
      }
    }
    edgesThroughCell.erase(id);
  }
  */
  return true;
}

template <class PointType, class DataType>
bool GridSearch< PointType, DataType>::Plan(GridPath< PointType, DataType >& path,
                                            bool removeDuplicateCells) {
  path.cells.clear();
  path.cost = 0;

  std::vector< CellEdge* > edgePath;

  if (Search< Cell<PointType, DataType> , GridPath<PointType, DataType>  >::Plan(edgePath) < 0)
    return false;

  typename vector< CellEdge* >::iterator it;
  for (auto&& edge : edgePath) {
    // insert in reverse if connectivity was expanded backwards
    if (edge->data.fwd)
      path.cells.insert(
          path.cells.end(), edge->data.cells.begin(), edge->data.cells.end());
    else
      path.cells.insert(
          path.cells.end(), edge->data.cells.rbegin(), edge->data.cells.rend());

    path.cost += edge->data.cost;
  }

  if (removeDuplicateCells && path.cells.size() > 1) {
    //    typename vector< TCell >::iterator cit;
    vector< TCell > cells;
    cells.push_back(path.cells.front());
    for (auto cit = path.cells.begin(); (cit + 1) != path.cells.end(); ++cit) {
      if ((cit->c - (cit + 1)->c).norm() > 1e-10)
        cells.push_back(*(cit + 1));
    }
    path.cells = cells;
  }
  // now add the very last one
  //    path.cells.push_back(edgePath.back()->data.cells.back());

  return true;
}


template < class PointType, class DataType >
void GridSearch< PointType, DataType >::OptPath(const GridPath< PointType, DataType >& path,
                                      GridPath< PointType, DataType >& optPath,
                                      double freeCost,
                                      double traceStep) const {
  double len = 0;

  optPath.cells.clear();
  optPath.cost = 0;

  if (path.cells.size() == 2) {
    optPath.cells = path.cells;
    optPath.cost = path.cost;
    return;
  }
  //  typename vector< Cell<PointType, DataType> >::const_iterator it0;
  auto it0 = path.cells.begin();
  //  typename vector< Cell<PointType, DataType> >::const_iterator it1;
  auto it1 = it0 + 1;

  PointType x0 = it0->c;
  PointType x1 = it1->c;

  PointType dx0 = x1 - x0;
  double dn = dx0.norm();
  dx0 /= dn;

  optPath.cells.push_back(path.cells[0]);

  if (traceStep <= 0)
    traceStep = grid.cs.norm();

  for (; it1 != path.cells.end() - 1; ++it1) {
    x1 = it1->c;
    PointType x2 = (it1 + 1)->c;
    PointType dx1 = x2 - x0;
    dn = dx1.norm();
    assert(dn > 1e-12);
    dx1 /= dn;

    if ((dx0 - dx1).norm() > 1e-16) {
      dn = (x0 - x2).norm();
      for (double d = traceStep; d < dn; d += traceStep) {
        PointType x = x0 + dx1 * d;
        int id = grid.Id(x);
        assert(id >= 0 && id < grid.nc);
        //        if (!grid.cells[id] || grid.cells[id]->cost > freeCost) {
        if (!grid.cells[id] || !connectivity.Free(grid.cells[id]->data)) {
          optPath.cells.push_back(*it1);
          x0 = x1;
          break;
        }
      }
      dx0 = x2 - x0;
      dn = dx0.norm();
      dx0 /= dn;
    }
  }

  optPath.cells.push_back(path.cells.back());
}

/*
template < int n, class CellData, class Tp >
void GridSearch< n, CellData, Tp >::SplinePath(const GridPath< n, CellData, Tp >& path,
                                         std::vector< PointType >& splinePath,
                                         double traceStep) const {
  std::vector< double > steps(path.cells.size());
  std::vector< std::vector< double > > pts(
      n, std::vector< double >(path.cells.size()));
  for (int i = 0; i < (int)(path.cells.size()); i++) {
    steps[i] = i;
    for (int j = 0; j < n; j++) {
      pts[j][i] = path.cells[i].c(j);
    }
  }

  std::vector< Spline< double, double > > sps;
  for (int i = 0; i < n; i++) {
    sps.push_back(Spline< double, double >(steps, pts[i]));
  }
  int count = (path.cells.size() - 1) / traceStep;
  splinePath.resize(count);
  // splineCells.cells.resize(count);
  // splineCells.cost = 0;

  for (int i = 0; i < count; i++) {
    PointType pti;
    for (int j = 0; j < n; j++) {
      pti(j) = sps[j][i * traceStep];
    }
    splinePath[i] = pti;

    // std::cout << i*traceStep << std::endl;
    // std::cout << pti.transpose() << std::endl;
    //  TODO(comment this out)
    Cell<n,CellData>* cell = grid.Get(splinePath[i]);
    if(!cell)
    {
      std::cout << "dsl::SplinePath: Cell does not exist" << std::endl;
      return;
    }
    splineCells.cells[i] = *(grid.Get(splinePath[i]));
    if(i > 0)
    {
      splineCells.cost +=
    (splineCells.cells[i].c-splineCells.cells[i-1].c).norm();
    }
  }
}
*/

template < class PointType, class DataType>
double GridSearch< PointType, DataType>::GetCost(const PointType& x) const {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return 0;

  const Cell< PointType, DataType >* cell = grid.Get(id);
  if (!cell) {
    std::cout << "[W] GridSearch::GetCost: no cell at position "
              << x.transpose() << std::endl;
    return 0;
  }
  return cell->cost;
}


template < class PointType, class DataType >
bool GridSearch< PointType, DataType>::SetCost(const PointType& x, double cost) {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return false;

  const Cell<PointType, DataType >* cell = grid.Get(id);
  if (!cell) {
    // std::cout << "[W] GridSearch::SetCost: no cell at position " <<
    // x.transpose() << std::endl;
    return false;
  }

  //    std::cout << "[W] GridSearch::SetCost: checking position " <<
  //    x.transpose() << std::endl;

  // if the cost has not changed simply return
  //  if (Search< Cell< PointType, DataType >, GridPath< PointType, DataType > >::Eq(cost, cell->cost))
  //    return false;

  //  cell->cost = cost;

  // change the costs of eall edges passing through this cell
  
  if (trackEdgesThroughCell) {
    typename map< int, vector< CellEdge* > >::iterator ceit =
        edgesThroughCell.find(id);
    if (ceit != edgesThroughCell.end()) {
      vector< CellEdge* > edges = ceit->second;
      //      std::cout << edges.size() << " edges passing through" << cell->c.transpose()
      //                << std::endl;

      typename vector< CellEdge* >::iterator ei;
      for (ei = edges.begin(); ei != edges.end(); ++ei) {
        this->ChangeCost(**ei, cost);
      }
    }
  }
 

  CellVertex* vertex = vertexMap[id];

  // if vertex hasn't been added yet, then we simply return successfully since
  // the search hasn't naturally explored this area previously
  if (!vertex)
    return true;

  // vertex->data.cost = cost;

  // fix edges are not tracked, then at best we can modify the once incoming and
  // outgoing from the cell
  // if edges are tracked this is already done above
  if (!trackEdgesThroughCell) {
    for (auto edge : vertex->in) {
      this->ChangeCost(*edge.second, cost);      
      //      this->ChangeCost(*ein->second, 
      //                       this->cost.Real(ein->second->from->data, vertex->data));
    }

    for (auto edge : vertex->out) {
      this->ChangeCost(*edge.second, cost);
      //      this->ChangeCost(*eout->second,
      //                       this->cost.Real(vertex->data, eout->second->to->data));
      
    }
  }
  return true;
}

template < class PointType, class DataType >
    bool GridSearch<PointType, DataType>::AddEdge(const PointType& x1, const PointType& x2) {
  CellVertex* from = GetVertex(x1);
  CellVertex* to = GetVertex(x2);
  if (!from || !to)
    return false;

  CellEdge* edge = new CellEdge(from, to, cost.Real(from->data, to->data));
  graph.AddEdge(*edge);
  return true;
}
}

#endif
