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
 *  1) GridSearch(grid, connectivity, cost)
 *  2) SetStart(Vector2d(x, y))
 *  3) SetGoal(Vector2d(x, y))
 *  4) Plan(path)
 *  5) follow path until map changes are observed
 *  6) for each change: SetCost(Vector2d(x, y), cost) 
 *  7) SetStart(Vector2d(x,y)) -- change the start to the current robot position
 *  8) goto 4
 *
 *
 *  Author: Marin Kobilarov
 */

namespace dsl {

using std::vector;

/**
 * Grid Planner that can compute the optimal path between two given cells on a
 *grid.
 * The path is a sequence of cells. Each cell in the grid has a cost and in
 *addition
 * there is cost of transitioning between two cells.
 *
 */
template < class PT, class DT, class ConnectionType = std::vector< Cell<PT, DT>* > >
  class GridSearch : public Search< Cell<PT, DT>, ConnectionType > {
public:
  using PointType = PT;
  using DataType = DT;

  using GridVertexData = Cell<PointType, DataType>;
  using GridEdgeData = ConnectionType;
  
  using CellVertex = Vertex< GridVertexData, GridEdgeData>;
  using CellEdge = Edge< GridVertexData, GridEdgeData>;

  using TypedCell = Cell<PointType, DataType>;
  using TypedCellPtr = shared_ptr<TypedCell>;

  using TypedGrid = Grid<PointType, DataType>;
    
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
  GridSearch(Grid< PointType, DataType >& grid,
             const GridConnectivity< PointType, DataType, ConnectionType >& connectivity,
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
  bool Plan(GridPath< PointType, DataType, ConnectionType >& path, bool removeDuplicateCells = true);

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

private:
  Graph< GridVertexData, GridEdgeData> graph; ///< the underlying graph

  Grid< PointType, DataType >& grid; ///< the grid
  const GridConnectivity< PointType, DataType, ConnectionType >&
      connectivity;              ///< the connectivity interface
  const GridCost< PointType, DataType >& cost; ///< the cost interface


  CellVertex** vertexMap; ///< vertex grid array

  bool trackEdgesThroughCell; ///< whether to keep a list of edges passing
  /// through each cell, this is useful for changing
  /// the costs these edges if the cost of that cell
  /// changes, or to remove these edges if the cell
  /// is removed
  std::map< int, std::vector< CellEdge* > >
      edgesThroughCell; ///< map of list of edges for each cell
};


template < class PointType, class DataType, class ConnectionType >
    GridSearch< PointType, DataType, ConnectionType>::GridSearch(
    Grid< PointType, DataType>& grid,
    const GridConnectivity< PointType, DataType, ConnectionType >& connectivity,
    const GridCost< PointType, DataType >& cost,
    bool expand,
    bool trackEdgesThroughCell)
    : Search< TypedCell, ConnectionType >(graph, cost),
    grid(grid),
    connectivity(connectivity),
    cost(cost),
    trackEdgesThroughCell(trackEdgesThroughCell) {

      vertexMap = new CellVertex* [grid.nc];
      memset(vertexMap, 0, grid.nc * sizeof(CellVertex*));

  if (expand) {
    for (int i = 0; i < grid.nc; ++i) {
      const TypedCellPtr cell = grid.Get(i);
      if (cell) {
        vertexMap[i] = new CellVertex(*cell);
        graph.AddVertex(*vertexMap[i]);
      }
    }

    // expand the successors of each vertex
    for (int i = 0; i < grid.nc; ++i) {
      const TypedCellPtr cell = grid.Get(i);
      if (cell) {
        CellVertex* from = vertexMap[i];
        assert(from);

        // generate successor paths
        std::vector< std::tuple<TypedCellPtr, ConnectionType, double> > paths;
        connectivity(*cell, paths);        

        for (auto&& path : paths) {
          // find cell where the end of the path falls
          TypedCellPtr toCell = std::get<0>(path);
          
          CellVertex* to = vertexMap[toCell->id];
          
          if (!to)
            continue;

          // path.cost = cost.Real(from->data, to->data);
          ConnectionType &connection = std::get<1>(path);
          double &cost = std::get<2>(path);
          
          CellEdge* edge = new CellEdge(connection, from, to, cost);
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

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType>::Expand(CellVertex& from, bool fwd) {
  // if this is true then all vertices have already been expanded at
  // construction
  if (fwd && from.succExpanded)
    return true;

  if (!fwd && from.predExpanded)
    return true;

  // cell must exist
  const TypedCellPtr cell = grid.Get(from.data.id);
  //  if (!cell) {
  //    std::cout << id << " " << from.data.c.transpose() << std::endl;    
  //  }
  assert(cell);
  
  std::vector< std::tuple<TypedCellPtr, ConnectionType, double> > paths;
  connectivity(*cell, paths, fwd);
  
  for (auto&& path : paths) {

    TypedCellPtr toCell = std::get<0>(path);
    assert(toCell);
    
    int id = toCell->id;

    // if this vertex doesn't exist, create it and add to graph
    if (!vertexMap[id]) {
      vertexMap[id] = new CellVertex(*toCell);
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

    //  path.cost = cost.Real(from.data, to->data);
    ConnectionType &connection = std::get<1>(path);
    double &cost = std::get<2>(path);
    
    CellEdge* edge = fwd ? new CellEdge(connection, &from, to, cost) :
                           new CellEdge(connection, to, &from, cost);
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

template < class PointType, class DataType, class ConnectionType >
    GridSearch< PointType, DataType, ConnectionType >::~GridSearch() {

  for (auto edge : graph.edges) {
    delete edge.second;
  }
  for (auto vertex : graph.vertices) {
    delete vertex.second;
  }

  delete[] vertexMap;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::SetStart(const PointType& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);

  // cell must exist
  const TypedCellPtr cell = grid.Get(id);
  if (!cell) {
    std::cout << "[W] GridSearch:SetStart: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }
  
  assert(cell->id == id);
  
  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  Search< Cell<PointType, DataType> , ConnectionType >::SetStart(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::SetGoal(const PointType& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetGoal: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);
  assert(id >= 0 && id < grid.nc);

  const TypedCellPtr cell = grid.Get(id);
  if (!cell) {
    std::cout << "[W] GridSearch:SetGoal: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }

  assert(cell->id == id);

  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  Search< Cell<PointType, DataType> , ConnectionType  >::SetGoal(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::RemoveCell(const PointType& x) {
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

  if (grid.cells[id]) {
    grid.cells[id].reset();
  } else {
    return false;
  }
  
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

template <class PointType, class DataType, class ConnectionType>
    bool GridSearch< PointType, DataType, ConnectionType>::Plan(GridPath< PointType, DataType, ConnectionType >& path,
                                                                bool removeDuplicateCells) {
  path.cells.clear();
  path.cost = 0;

  std::vector< CellEdge* > edgePath;

  if (Search< TypedCell, ConnectionType >::Plan(edgePath) < 0)
    return false;

  for (auto&& edge : edgePath) {
    // insert in reverse if connectivity was expanded backwards

    /* TODO: marin (this should not be necessary, unless the ConnectionType's are also reversed
    if (edge->data.fwd)
      path.cells.insert(
          path.cells.end(), edge->data.cells.begin(), edge->data.cells.end());
    else
      path.cells.insert(
          path.cells.end(), edge->data.cells.rbegin(), edge->data.cells.rend());
    */
    path.cells.push_back(edge->from->data);
    path.connections.push_back(edge->data);
    path.cost += edge->cost;
  }
  path.cells.push_back(edgePath.back()->to->data);


  if (removeDuplicateCells && path.cells.size() > 1) {
    std::vector< TypedCell > cells;
    cells.push_back(path.cells.front());
    for (auto cit = path.cells.begin(); (cit + 1) != path.cells.end(); ++cit) {
      if ((cit->c - (cit + 1)->c).norm() > 1e-10)
        cells.push_back(*(cit + 1));
    }
    path.cells = cells;
  }
  // now add the very last one: as per above not necessary any more
  //    path.cells.push_back(edgePath.back()->data.cells.back());

  return true;
}


template < class PointType, class DataType, class ConnectionType>
    double GridSearch< PointType, DataType, ConnectionType>::GetCost(const PointType& x) const {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return 0;

  Cell< PointType, DataType >* cell = grid.Get(id);
  if (!cell) {
    std::cout << "[W] GridSearch::GetCost: no cell at position "
              << x.transpose() << std::endl;
    return 0;
  }
  return cell->cost;
}


template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType>::SetCost(const PointType& x, double cost) {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return false;

  TypedCellPtr cell = grid.Get(id);
  if (!cell) {
    // std::cout << "[W] GridSearch::SetCost: no cell at position " <<
    // x.transpose() << std::endl;
    return false;
  }
  assert(id == cell->id);

  // change the costs of eall edges passing through this cell  
  if (trackEdgesThroughCell) {
    auto ceit = edgesThroughCell.find(id);
    if (ceit != edgesThroughCell.end()) {
      std::vector< CellEdge* > edges = ceit->second;
      for (auto&& edge : edges) {
        this->ChangeCost(*edge, cost);
      }
    }
  }
 

  CellVertex* vertex = vertexMap[id];

  // if vertex hasn't been added yet, then we simply return successfully since
  // the search hasn't naturally explored this area previously
  if (!vertex)
    return true;

  // fix edges are not tracked, then at best we can modify the once incoming and
  // outgoing from the cell
  // if edges are tracked this is already done above
  if (!trackEdgesThroughCell) {
    for (auto edge : vertex->in) {
      this->ChangeCost(*edge.second, cost);      
    }

    for (auto edge : vertex->out) {
      this->ChangeCost(*edge.second, cost);
    }
  }
  return true;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch<PointType, DataType, ConnectionType>::AddEdge(const PointType& x1, const PointType& x2) {
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
