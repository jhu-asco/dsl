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
template < class PointT,
           class DataT,
           class ConnectionT = std::vector< Cell< PointT, DataT >* > >
class GridSearch {
public:
  using VertexDataT = Cell< PointT, DataT >;
  using EdgeDataT = ConnectionT;

  using VertexT = Vertex< VertexDataT, EdgeDataT >;
  using EdgeT = Edge< VertexDataT, EdgeDataT >;

  using CellT = Cell< PointT, DataT >;
  using GridT = Grid< PointT, DataT >;

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
   * @param track_edges_through_cell whether to keep a list of edges passing
   * through each cell, this is useful for changing the costs these edges if the
   * cost of that cell changes, or to remove these edges if the cell is removed
   */
  GridSearch(const GridT& grid,
             const GridConnectivity< PointT, DataT, ConnectionT >& connectivity,
             const GridCost< PointT, DataT >& cost,
             Method method = Method::kDstar,
             bool expand = false,
             bool track_edges_through_cell = false);

  virtual ~GridSearch();

  /**
   * Change the cost of an individual cell
   * @param x position
   * @param cost cost
   */
  bool setCost(const PointT& x, double cost);

  /**
   * Get the cost of an individual cell
   * @param x position
   * @return cost
   */
  double getCost(const PointT& x) const;

  /**
   * Set start location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool setStart(const PointT& x);

  /**
   * Set goal location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool addGoal(const PointT& x);

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
  bool expand(VertexT& from, bool fwd = true);

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
  bool plan(GridPath< PointT, DataT, ConnectionT >& path,
            bool remove_duplicate_cells = true);

  /**
   * Useful method to get the graph vertex at position x
   * @param x position
   * @return corresponding vertex or 0 if none there
   */
  VertexT* getVertex(const PointT& x) const;

  /**
   * Useful method to remove a vertex at position x
   * @param x Euclidean point vector
   */
  bool removeCell(const PointT& x);

  /**
   * Useful method for adding edges b/n vertices
   * @param x1 start point
   * @param x2 end point
   */
  bool addEdge(const PointT& x1, const PointT& x2);

  Search< VertexDataT, EdgeDataT >* search = 0;

  const Graph< VertexDataT, EdgeDataT >& getGraph() {
    return graph;
  }

private:
  Graph< VertexDataT, EdgeDataT > graph; ///< the underlying graph

  const GridT& grid; ///< the grid
  const GridConnectivity< PointT, DataT, ConnectionT >&
      connectivity;                      ///< the connectivity interface
  const GridCost< PointT, DataT >& cost; ///< the cost interface

  VertexT** vertex_map; ///< vertex grid array

  bool track_edges_through_cell; ///< whether to keep a list of edges passing
  /// through each cell, this is useful for changing
  /// the costs these edges if the cost of that cell
  /// changes, or to remove these edges if the cell
  /// is removed
  std::map< int, std::vector< EdgeT* > >
      edges_through_cell; ///< map of list of edges for each cell
};

template < class PointT, class DataT, class ConnectionT >
GridSearch< PointT, DataT, ConnectionT >::GridSearch(
    const Grid< PointT, DataT >& grid,
    const GridConnectivity< PointT, DataT, ConnectionT >& connectivity,
    const GridCost< PointT, DataT >& cost,
    Method method,
    bool expand,
    bool track_edges_through_cell)
  : grid(grid),
    connectivity(connectivity),
    cost(cost),
    track_edges_through_cell(track_edges_through_cell) {
  switch(method) {
    case Method::kDstar:
      search = new Dstar< VertexDataT, EdgeDataT >(graph, cost);
      break;
    case Method::kLpAstar:
      search = new LpAstar< VertexDataT, EdgeDataT >(graph, cost);
      break;
    default:
      std::cerr << "[E] GridSearch::GridSearch: wrong type" << (int)method;
      return;
  }

  auto expand_callback =
      [this](VertexT& from, bool fwd) { return this->expand(from, fwd); };
  search->setExpandCallback(expand_callback);

  vertex_map = new VertexT* [grid.nc];
  memset(vertex_map, 0, grid.nc * sizeof(VertexT*));

  if (expand) {
    for (int i = 0; i < grid.nc; ++i) {
      const CellT* cell = grid.data(i);
      if (cell) {
        vertex_map[i] = new VertexT(*cell);
        graph.addVertex(*vertex_map[i]);
      }
    }

    // expand the successors of each vertex
    for (int i = 0; i < grid.nc; ++i) {
      const CellT* cell = grid.data(i);
      if (cell) {
        VertexT* from = vertex_map[i];
        assert(from);

        // generate successor paths
        std::vector< std::tuple< CellT*, ConnectionT, double > > paths;
        connectivity(*cell, paths);

        for (auto&& path : paths) {
          // find cell where the end of the path falls
          CellT* to_cell = std::get< 0 >(path);

          VertexT* to = vertex_map[to_cell->id];

          if (!to)
            continue;

          // path.cost = cost.real(from->data, to->data);
          ConnectionT& connection = std::get< 1 >(path);
          double &cost = std::get<2>(path);

          EdgeT* edge = new EdgeT(connection, from, to, cost);
          graph.addEdge(*edge);

          /*
          TODO(marin): re-enable this functionality if it is needed
          if (track_edges_through_cell) {
            // iterate through path and add the new edge to the list
            // of edges passing through each cell of the path
            for (auto&& cells : path.cells) {
              int id = grid.computeId(cells.first);
              //              typename map< int, vector< EdgeT* > >::iterator
          etcit;
              auto etcit = edges_through_cell.find(id);
              if (etcit == edges_through_cell.end()) {
                vector< EdgeT* > edges;
                edges.push_back(edge);
                edges_through_cell[id] = edges;
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

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::expand(VertexT& from, bool fwd) {
  // if this is true then all vertices have already been expanded at
  // construction
  if (fwd && from.succ_expanded)
    return true;

  if (!fwd && from.pred_expanded)
    return true;

  // cell must exist
  const Cell< PointT, DataT >* cell = grid.data(from.data.id);
  assert(cell);

  std::vector< std::tuple< CellT*, ConnectionT, double > > paths;
  connectivity(*cell, paths, fwd);

  for (auto&& path : paths) {
    CellT* to_cell = std::get< 0 >(path);
    assert(to_cell);

    int id = to_cell->id;

    // if this vertex doesn't exist, create it and add to graph
    if (!vertex_map[id]) {
      vertex_map[id] = new VertexT(*to_cell);
      graph.addVertex(*vertex_map[id]);
    }
    VertexT* to = vertex_map[id];

    // if fwd and incoming edge from->to exists, then do not create a new one
    if (from.find(*to, !fwd))
      continue;

    ConnectionT& connection = std::get< 1 >(path);
    double &cost = std::get<2>(path);

    EdgeT* edge = fwd ? new EdgeT(connection, &from, to, cost) :
                        new EdgeT(connection, to, &from, cost);
    graph.addEdge(*edge);

    /*
    TODO(marin): re-enable this functionality if it is needed
    if (track_edges_through_cell) {
      // iterate through path and add the new edge to the list
      // of edges passing through each cell of the path
      for (auto cell : path.cell) {
        int id = grid.computeId(cell.first);
        typename map< int, vector< EdgeT* > >::iterator etcit;
        etcit = edges_through_cell.find(id);
        if (etcit == edges_through_cell.end()) {
          vector< EdgeT* > edges;
          edges.push_back(edge);
          edges_through_cell[id] = edges;
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

template < class PointT, class DataT, class ConnectionT >
GridSearch< PointT, DataT, ConnectionT >::~GridSearch() {
  for (auto edge : graph.edges) {
    delete edge.second;
  }
  for (auto vertex : graph.vertices) {
    delete vertex.second;
  }

  delete[] vertex_map;
  delete search;
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::setStart(const PointT& x) {
  if (!grid.valid(x)) {
    std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.computeId(x);

  // cell must exist
  const Cell< PointT, DataT >* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch:SetStart: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }

  assert(cell->id == id);

  // if it's not added previously add it
  if (!vertex_map[id]) {
    vertex_map[id] = new VertexT(*cell);
    graph.addVertex(*vertex_map[id]);
  }

  search->setStart(*vertex_map[id]);

  return true;
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::addGoal(const PointT& x) {
  if (!grid.valid(x)) {
    std::cout << "[W] GridSearch:setGoal: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.computeId(x);
  assert(id >= 0 && id < grid.nc);

  const CellT* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch:setGoal: cell at=" << x.transpose() << " does not exist!"
              << std::endl;
    return false;
  }

  assert(cell->id == id);

  // if it's not added previously add it
  if (!vertex_map[id]) {
    vertex_map[id] = new VertexT(*cell);
    graph.addVertex(*vertex_map[id]);
  }

  search->addGoal(*vertex_map[id]);

  return true;
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::removeCell(const PointT& x) {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    std::cout << "[W] GridSearch::removeCell: id=" << id << " out of bounds!" << std::endl;
    return false;
  }

  VertexT* v = vertex_map[id];
  if (v) {
    graph.removeVertex(*v);
    delete v;
    vertex_map[id] = 0;
  } else {
    std::cout << "[W] GridSearch::removeCell: vertex with id=" << id << " doesnt exist!" << std::endl;
    return false;
  }

  if (grid.values[id]) {
    delete grid.values[id];
    grid.values[id] = 0;
  } else {
    return false;
  }

  /*
  TODO(marin): re-enable this functionality if it is needed
  // if edges through cells are tracked, then go through them and remove them
  if (track_edges_through_cell) {
    typename map< int, vector< EdgeT* > >::iterator ceit =
        edges_through_cell.find(id);
    if (ceit != edges_through_cell.end()) {
      vector< EdgeT* > edges = ceit->second;
      typename vector< EdgeT* >::iterator ei;
      for (ei = edges.begin(); ei != edges.end(); ++ei) {
        this->graph.removeEdge(**ei);
      }
    }
    edges_through_cell.erase(id);
  }
  */
  return true;
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::plan(
    GridPath< PointT, DataT, ConnectionT >& path, bool remove_duplicate_cells) {
  path.cells.clear();
  path.cost = 0;

  std::vector< EdgeT* > edgePath;

  if (search->plan(edgePath) < 0)
    return false;

  for (auto&& edge : edgePath) {
    // insert in reverse if connectivity was expanded backwards

    /* TODO(marin) (this should not be necessary, unless the ConnectionT's are
    also reversed (e.g. if the primitives are generated backwards for some
    reason)
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
    std::vector< CellT > cells;
    cells.push_back(path.cells.front());
    for (auto cit = path.cells.begin(); (cit + 1) != path.cells.end(); ++cit) {
      if ((cit->center - (cit + 1)->center).norm() > 1e-10)
        cells.push_back(*(cit + 1));
    }
    path.cells = cells;
  }
  // now add the very last one: as per above not necessary any more
  //    path.cells.push_back(edgePath.back()->data.cells.back());

  return true;
}

template < class PointT, class DataT, class ConnectionT >
double
    GridSearch< PointT, DataT, ConnectionT >::getCost(const PointT& x) const {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return 0;
  }

  Cell< PointT, DataT >* cell = grid.data(id);
  if (!cell) {
    std::cout << "[W] GridSearch::getCost: no cell at position "
              << x.transpose() << std::endl;
    return 0;
  }
  return cell->cost;
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::setCost(const PointT& x,
                                                       double cost) {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return false;
  }

  Cell< PointT, DataT >* cell = grid.data(id);
  if (!cell) {
    return false;
  }
  assert(id == cell->id);

  // change the costs of eall edges passing through this cell
  if (track_edges_through_cell) {
    auto ceit = edges_through_cell.find(id);
    if (ceit != edges_through_cell.end()) {
      std::vector< EdgeT* > edges = ceit->second;
      for (auto& edge : edges) {
        search->changeCost(*edge, cost);
      }
    }
  }

  VertexT* vertex = vertex_map[id];

  // if vertex hasn't been added yet, then we simply return successfully since
  // the search hasn't naturally explored this area previously
  if (!vertex) {
    return true;
  }

  // fix edges are not tracked, then at best we can modify the once incoming and
  // outgoing from the cell
  // if edges are tracked this is already done above
  if (!track_edges_through_cell) {
    for (auto& edge : vertex->in) {
      search->changeCost(*edge.second, cost);
    }

    for (auto& edge : vertex->out) {
      search->changeCost(*edge.second, cost);
    }
  }
  return true;
}

template < class PointT, class DataT, class ConnectionT >
Vertex< Cell< PointT, DataT >, ConnectionT >*
    GridSearch< PointT, DataT, ConnectionT >::getVertex(const PointT& x) const {
  int id = grid.computeId(x);
  if (id < 0 || id >= grid.nc) {
    return 0;
  }
  return vertex_map[id];
}

template < class PointT, class DataT, class ConnectionT >
bool GridSearch< PointT, DataT, ConnectionT >::addEdge(const PointT& x1,
                                                       const PointT& x2) {
  VertexT* from = getVertex(x1);
  VertexT* to = getVertex(x2);
  if (!from || !to)
    return false;

  EdgeT* edge = new EdgeT(from, to, cost.real(from->data, to->data));
  graph.addEdge(*edge);
  return true;
}
}
