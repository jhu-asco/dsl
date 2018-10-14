// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <iostream>
#include "lpastar.h"
#include "dstar.h"
#include "gridcost.h"
#include "gridconnectivity.h"
#include "gridpath.h"
#include "grid.h"

/**
 *  Grid based D*-Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base LpAstar class
 *
 *  Supports n-dimensional grid with abstract "connectivity" interface
 *
 *  For example in 2D, the cell transition costs can be encoded as the euclidean
 *distance
 *  b/n the centrs of the cells (i.e. each cell has 8 neighbors:
 *  the neighbors at N,S,W,E have transition cost of 1, and the
 *  neighbors at NE, SE, NW, SW have transition costs of sqrt(2)
 *  These transition costs are added to the maximum of the values of
 *  two neighboring cells to compute the cost of their connecting edge.
 *
 *
 *  The planner is used as follows:
 *
 *  1) GridSearch(grid, connectivity, cost)
 *  2) setStart(Vector2d(x, y))
 *  3) setGoal(Vector2d(x, y))
 *  4) plan(path)
 *  5) follow path until map changes are observed
 *  6) for each change: setCost(Vector2d(x, y), cost)
 *  7) setStart(Vector2d(x,y)) -- change the start to the current robot position
 *  8) goto 4
 *
 *
 *  Author: Marin Kobilarov
 */

namespace dsl {

/**
 * Grid Planner that can compute the optimal path between two given cells on a
 *grid.
 * The path is a sequence of cells. Each cell in the grid has a cost and in
 *addition
 * there is cost of transitioning between two cells.
 *
 */
template < class PointType, class DataType, class ConnectionType = std::vector< Cell<PointType, DataType>* > >
  class GridSearch  {
public:

  using GridVertexData = Cell<PointType, DataType>;
  using GridEdgeData = ConnectionType;

  using CellVertex = Vertex< GridVertexData, GridEdgeData>;
  using CellEdge = Edge< GridVertexData, GridEdgeData>;

  using CellType = Cell<PointType, DataType>;
  using GridType = Grid<PointType, DataType>;

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
   * @param method method (Dstar or LpAstar)
   * @param expand whether to construct/expand the whole graph at init
   * @param trackEdgesThroughCell whether to keep a list of edges passing
   * through each cell, this is useful for changing the costs these edges if the
   * cost of that cell changes, or to remove these edges if the cell is removed
   */
  GridSearch(const Grid< PointType, DataType >& grid,
             const GridConnectivity< PointType, DataType, ConnectionType >& connectivity,
             const GridCost< PointType, DataType >& cost,
             Method method = Method::kDstar,
             bool expand = false,
             bool trackEdgesThroughCell = false);

  virtual ~GridSearch();

  /**
   * Change the cost of an individual cell
   * @param x position
   * @param cost cost
   */
  bool setCost(const PointType& x, double cost);

  /**
   * Get the cost of an individual cell
   * @param x position
   * @return cost
   */
  double getCost(const PointType& x) const;

  /**
   * Set start location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool setStart(const PointType& x);

  /**
   * Set goal location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool addGoal(const PointType& x);

  /**
   * Expand successors or predecessors of a given vertex. This is mostly used
   * internally
   * during search-> In some cases can be "pre-called" externally in advance of
   * searching in this area of
   * the grid, for efficiency purposes.
   * @param from the given vertex
   * @param fwd if true then expand forward in time, i.e. successors, otherwise
   * expand predecessors
   * @return true on success
   */
  bool expand(CellVertex& from, bool fwd = true);

  /**
   * Compute path b/n start and goal vertices
   * these vertices should be already set
   * @param path the resulting path
   * @param remove_duplicate_cells remove duplicate consecutive cells (this might
   * appear if two consecutive edges
   * contain cells in a way that the last cell of the first edge overlaps with
   * the first cell of the next edge)
   * @return true on success
   */
  bool plan(GridPath< PointType, DataType, ConnectionType >& path, bool remove_duplicate_cells = true);

  /**
   * Useful method to get the graph vertex at position x
   * @param x position
   * @return corresponding vertex or 0 if none there
   */
  CellVertex* getVertex(const PointType& x) const;

  /**
   * Useful method to remove a vertex at position x
   * @param x Euclidean point vector
   */
  bool removeCell(const PointType& x);

  /**
   * Useful method for adding edges b/n vertices
   * @param x1 start point
   * @param x2 end point
   */
  bool addEdge(const PointType& x1, const PointType& x2);

  Search<GridVertexData, GridEdgeData> *search;

  const Graph< GridVertexData, GridEdgeData>& GetGraph() {
    return graph;
  }

private:
  Graph< GridVertexData, GridEdgeData> graph; ///< the underlying graph


  const Grid< PointType, DataType >& grid; ///< the grid
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
    const Grid< PointType, DataType>& grid,
    const GridConnectivity< PointType, DataType, ConnectionType >& connectivity,
    const GridCost< PointType, DataType >& cost,
    Method method,
    bool expand,
    bool trackEdgesThroughCell)
    : grid(grid),
    connectivity(connectivity),
    cost(cost),
    trackEdgesThroughCell(trackEdgesThroughCell) {

  switch(method) {
    case Method::kDstar:
      search = new Dstar<GridVertexData, GridEdgeData>(graph, cost);
      break;
    case Method::kLpAstar:
      search = new LpAstar<GridVertexData, GridEdgeData>(graph, cost);
      break;
    default:
      std::cerr << "[E] GridSearch::GridSearch: wront type" << (int)method;
      return;
  }

  auto expand_callback =
    [this](CellVertex& from, bool fwd) {return this->expand(from, fwd);};
  search->setExpandCallback(expand_callback);

  vertexMap = new CellVertex* [grid.nc];
  memset(vertexMap, 0, grid.nc * sizeof(CellVertex*));

  if (expand) {
    for (int i = 0; i < grid.nc; ++i) {
      const CellType *cell = grid.data(i);
      if (cell) {
        vertexMap[i] = new CellVertex(*cell);
        graph.addVertex(*vertexMap[i]);
      }
    }

    // expand the successors of each vertex
    for (int i = 0; i < grid.nc; ++i) {
      const CellType *cell = grid.data(i);
      if (cell) {
        CellVertex* from = vertexMap[i];
        assert(from);

        // generate successor paths
        std::vector< std::tuple<CellType*, ConnectionType, double> > paths;
        connectivity(*cell, paths);

        for (auto&& path : paths) {
          // find cell where the end of the path falls
          CellType *to_cell = std::get<0>(path);

          CellVertex* to = vertexMap[to_cell->id];

          if (!to)
            continue;

          // path.cost = cost.real(from->data, to->data);
          ConnectionType &connection = std::get<1>(path);
          double &cost = std::get<2>(path);

          CellEdge* edge = new CellEdge(connection, from, to, cost);
          graph.addEdge(*edge);

          /*
          TODO(marin): re-enable this functionality if it is needed
          if (trackEdgesThroughCell) {
            // iterate through path and add the new edge to the list
            // of edges passing through each cell of the path
            for (auto&& cells : path.cells) {
              int id = grid.computeId(cells.first);
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
        from->pred_expanded = true;
        from->succ_expanded = true;
      }
    }
  }
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType>::expand(CellVertex& from, bool fwd) {
  // if this is true then all vertices have already been expanded at
  // construction
  if (fwd && from.succ_expanded)
    return true;

  if (!fwd && from.pred_expanded)
    return true;

  // cell must exist
  const Cell<PointType, DataType>* cell = grid.data(from.data.id);
  assert(cell);

  std::vector< std::tuple<CellType*, ConnectionType, double> > paths;
  connectivity(*cell, paths, fwd);

  for (auto&& path : paths) {

    CellType *to_cell = std::get<0>(path);
    assert(to_cell);

    int id = to_cell->id;

    // if this vertex doesn't exist, create it and add to graph
    if (!vertexMap[id]) {
      vertexMap[id] = new CellVertex(*to_cell);
      graph.addVertex(*vertexMap[id]);
    }
    CellVertex* to = vertexMap[id];

    // if fwd and incoming edge from->to exists, then do not create a new one
    if (from.find(*to, !fwd))
      continue;

    ConnectionType &connection = std::get<1>(path);
    double &cost = std::get<2>(path);

    CellEdge* edge = fwd ? new CellEdge(connection, &from, to, cost) :
                           new CellEdge(connection, to, &from, cost);
    graph.addEdge(*edge);

    /*
    TODO(marin): re-enable this functionality if it is needed
    if (trackEdgesThroughCell) {
      // iterate through path and add the new edge to the list
      // of edges passing through each cell of the path
      for (auto cell : path.cell) {
        int id = grid.computeId(cell.first);
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
    from.succ_expanded = true;
  else
    from.pred_expanded = true;

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
  delete search;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::setStart(const PointType& x) {
  if (!grid.valid(x)) {
    std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.computeId(x);

  // cell must exist
  const Cell<PointType, DataType >* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch:SetStart: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }

  assert(cell->id == id);

  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.addVertex(*vertexMap[id]);
  }

  search->setStart(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::addGoal(const PointType& x) {
  if (!grid.valid(x)) {
    std::cout << "[W] GridSearch:setGoal: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.computeId(x);
  assert(id >= 0 && id < grid.nc);

  const CellType* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch:setGoal: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }

  assert(cell->id == id);

  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.addVertex(*vertexMap[id]);
  }

  search->addGoal(*vertexMap[id]);

  return true;
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType >::removeCell(const PointType& x) {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    std::cout << "[W] GridSearch::removeCell: id=" << id << " out of bounds!" << std::endl;
    return false;
  }

  CellVertex* v = vertexMap[id];
  if (v) {
    graph.removeVertex(*v);
    delete v;
    vertexMap[id] = 0;
  } else {
    std::cout << "[W] GridSearch::removeCell: vertex with id=" << id << " doesnt exist!" << std::endl;
    return false;
  }

  if (grid.cells[id]) {
    delete grid.cells[id];
    grid.cells[id] = 0;
  } else {
    return false;
  }

  /*
  TODO(marin): re-enable this functionality if it is needed
  // if edges through cells are tracked, then go through them and remove them
  if (trackEdgesThroughCell) {
    typename map< int, vector< CellEdge* > >::iterator ceit =
        edgesThroughCell.find(id);
    if (ceit != edgesThroughCell.end()) {
      vector< CellEdge* > edges = ceit->second;
      typename vector< CellEdge* >::iterator ei;
      for (ei = edges.begin(); ei != edges.end(); ++ei) {
        this->graph.removeEdge(**ei);
      }
    }
    edgesThroughCell.erase(id);
  }
  */
  return true;
}

template <class PointType, class DataType, class ConnectionType>
    bool GridSearch< PointType, DataType, ConnectionType>::plan(GridPath< PointType, DataType, ConnectionType >& path,
                                                                bool remove_duplicate_cells) {
  path.cells.clear();
  path.cost = 0;

  std::vector< CellEdge* > edgePath;

  if (search->plan(edgePath) < 0)
    return false;

  for (auto&& edge : edgePath) {
    // insert in reverse if connectivity was expanded backwards

    /* TODO(marin) (this should not be necessary, unless the ConnectionType's are also reversed (e.g. if the primitives are generated backwards for some reason)
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


  if (remove_duplicate_cells && path.cells.size() > 1) {
    std::vector< CellType > cells;
    cells.push_back(path.cells.front());
    for (auto cit = path.cells.begin(); (cit + 1) != path.cells.end(); ++cit) {
      if ((cit->centr - (cit + 1)->centr).norm() > 1e-10)
        cells.push_back(*(cit + 1));
    }
    path.cells = cells;
  }
  // now add the very last one: as per above not necessary any more
  //    path.cells.push_back(edgePath.back()->data.cells.back());

  return true;
}


template < class PointType, class DataType, class ConnectionType>
    double GridSearch< PointType, DataType, ConnectionType>::getCost(const PointType& x) const {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return 0;
  }

  Cell< PointType, DataType >* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch::getCost: no cell at position "
              << x.transpose() << std::endl;
    return 0;
  }
  return cell->cost;
}


template < class PointType, class DataType, class ConnectionType >
    bool GridSearch< PointType, DataType, ConnectionType>::setCost(const PointType& x, double cost) {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return false;
  }

  Cell<PointType, DataType >* cell = grid.data(id);
  if (!cell) {
    return false;
  }
  assert(id == cell->id);

  // change the costs of eall edges passing through this cell
  if (trackEdgesThroughCell) {
    auto ceit = edgesThroughCell.find(id);
    if (ceit != edgesThroughCell.end()) {
      std::vector< CellEdge* > edges = ceit->second;
      for (auto& edge : edges) {
        search->changeCost(*edge, cost);
      }
    }
  }


  CellVertex* vertex = vertexMap[id];

  // if vertex hasn't been added yet, then we simply return successfully since
  // the search hasn't naturally explored this area previously
  if (!vertex) {
    return true;
  }

  // fix edges are not tracked, then at best we can modify the once incoming and
  // outgoing from the cell
  // if edges are tracked this is already done above
  if (!trackEdgesThroughCell) {
    for (auto& edge : vertex->in) {
      search->changeCost(*edge.second, cost);
    }

    for (auto& edge : vertex->out) {
      search->changeCost(*edge.second, cost);
    }
  }
  return true;
}

template < class PointType, class DataType, class ConnectionType >
    Vertex< Cell<PointType, DataType>, ConnectionType>* GridSearch<PointType, DataType, ConnectionType>::getVertex(const PointType& x) const {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return 0;
  }
  return vertexMap[id];
}

template < class PointType, class DataType, class ConnectionType >
    bool GridSearch<PointType, DataType, ConnectionType>::addEdge(const PointType& x1, const PointType& x2) {
  CellVertex* from = getVertex(x1);
  CellVertex* to = getVertex(x2);
  if (!from || !to)
    return false;

  CellEdge* edge = new CellEdge(from, to, cost.real(from->data, to->data));
  graph.addEdge(*edge);
  return true;
}
}
