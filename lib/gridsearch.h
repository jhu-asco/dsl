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

/**
 * Grid Planner that can compute the optimal path between two given cells on a
 *grid.
 * The path is a sequence of cells. Each cell in the grid has a cost and in
 *addition
 * there is cost of transitioning between two cells.
 *
 */
template < int n, class Tc = Matrix< double, n, 1 >, class Tp = vector< Tc > >
class GridSearch : public Search< Cell< n, Tc >, GridPath< n, Tc, Tp > > {
public:
  typedef Matrix< double, n, 1 > Vectornd;
  typedef Matrix< int, n, 1 > Vectorni;
  typedef Vertex< Cell< n, Tc >, GridPath< n, Tc, Tp > > CellVertex;
  typedef Edge< Cell< n, Tc >, GridPath< n, Tc, Tp > > CellEdge;

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
  GridSearch(const Grid< n, Tc >& grid,
             const GridConnectivity< n, Tc, Tp >& connectivity,
             const GridCost< n, Tc >& cost,
             bool expand = false,
             bool trackEdgesThroughCell = false);

  virtual ~GridSearch();

  /**
   * Change the cost of an individual cell
   * @param x position
   * @param cost cost
   */
  bool SetCost(const Vectornd& x, double cost);

  /**
   * Get the cost of an individual cell
   * @param x position
   * @return cost
   */
  double GetCost(const Vectornd& x) const;

  //    void SetMap(const double *map);

  /**
   * Set start location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool SetStart(const Vectornd& x);

  /**
   * Set goal location in the map
   * @param x euclidean point vector
   * @return true on success
   */
  bool SetGoal(const Vectornd& x);

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
  bool Plan(GridPath< n, Tc, Tp >& path, bool removeDuplicateCells = true);

  /**
   * Useful method to get the graph vertex at position x
   * @param x position
   * @return corresponding vertex or 0 if none there
   */
  CellVertex* GetVertex(const Vectornd& x) const;

  /**
   * Useful method to remove a vertex at position x
   * @param x Euclidean point vector
   */
  bool RemoveCell(const Vectornd& x);

  /**
   * Useful method for adding edges b/n vertices
   * @param x1 start point
   * @param x2 end point
   */
  bool AddEdge(const Vectornd& x1, const Vectornd& x2);

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
  void OptPath(const GridPath< n, Tc, Tp >& path,
               GridPath< n, Tc, Tp >& optPath,
               double freeCost = 1e-3,
               double traceStep = -1.0) const;

  /**
   * Experimental path "smoothing" function
   * @param path original path
   * @param splinePath spline path
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is 0.1)
   */
  void SplinePath(const GridPath< n, Tc, Tp >& path,
                  std::vector< Vectornd >& splinePath,
                  // GridPath<n,Tc> &splineCells,
                  double traceStep = 0.1) const;

protected:
  Graph< Cell< n, Tc >, GridPath< n, Tc, Tp > > graph; ///< the underlying graph

  const Grid< n, Tc >& grid; ///< the grid
  const GridConnectivity< n, Tc, Tp >&
      connectivity;              ///< the connectivity interface
  const GridCost< n, Tc >& cost; ///< the cost interface

  CellVertex** vertexMap; ///< vertex grid array

  bool trackEdgesThroughCell; ///< whether to keep a list of edges passing
  /// through each cell, this is useful for changing
  /// the costs these edges if the cost of that cell
  /// changes, or to remove these edges if the cell
  /// is removed
  std::map< int, vector< CellEdge* > >
      edgesThroughCell; ///< map of list of edges for each cell
};

template < int n, class Tc, class Tp >
GridSearch< n, Tc, Tp >::GridSearch(
    const Grid< n, Tc >& grid,
    const GridConnectivity< n, Tc, Tp >& connectivity,
    const GridCost< n, Tc >& cost,
    bool expand,
    bool trackEdgesThroughCell)
  : Search< Cell< n, Tc >, GridPath< n, Tc, Tp > >(graph, cost),
    grid(grid),
    connectivity(connectivity),
    cost(cost),
    trackEdgesThroughCell(trackEdgesThroughCell) {
  vertexMap = new CellVertex* [grid.nc];
  memset(vertexMap, 0, grid.nc * sizeof(CellVertex*));

  if (expand) {
    // initialize and add all vertices
    for (int i = 0; i < grid.nc; ++i) {
      if (grid.cells[i]) {
        vertexMap[i] = new CellVertex(*grid.cells[i]);
        graph.AddVertex(*vertexMap[i]);
      }
    }

    // expand the successors of each vertex
    for (int i = 0; i < grid.nc; ++i) {
      if (grid.cells[i]) {
        CellVertex* from = vertexMap[i];
        assert(from);

        // generate successor paths
        std::vector< GridPath< n, Tc, Tp > > paths;
        connectivity(*grid.cells[i], paths);

        typename std::vector< GridPath< n, Tc, Tp > >::iterator it;
        for (it = paths.begin(); it != paths.end(); ++it) {
          GridPath< n, Tc, Tp >& path = *it;

          // find cell where the end of the path falls
          int id = grid.Id(path.cells.back().c);
          assert(id >= 0 && id < grid.nc);
          CellVertex* to = vertexMap[id];
          if (!to)
            continue;

          /*
          if (path.cells.size() >= 2) {
            path.cost = 0;
            typename vector<Cell<n, Tc> >::iterator cit = path.cells.begin();
            for (; (cit+1) != path.cells.end(); ++cit) {
              path.cost += cost.Real(*cit, *(cit+1));
            }
          } else {
            path.cost = cost.Real(from->data, to->data);
          }
          */

          //          path.cost = cost.Real(from->data, to->data);

          CellEdge* edge = new CellEdge(path, from, to, path.cost);
          graph.AddEdge(*edge);

          if (trackEdgesThroughCell) {
            // iterate through path and add the new edge to the list
            // of edges passing through each cell of the path
            for (typename vector< Cell< n, Tc > >::iterator cit =
                     path.cells.begin();
                 cit != path.cells.end();
                 ++cit) {
              int id = grid.Id(cit->c);
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
        }
        from->predExpanded = true;
        from->succExpanded = true;
      }
    }
  }
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::Expand(CellVertex& from, bool fwd) {
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
  Cell< n, Tc >* cell = grid.cells[id];
  assert(cell);

  std::vector< GridPath< n, Tc, Tp > > paths;
  connectivity(*cell, paths, fwd);

  for (int j = 0; j < paths.size(); ++j) {
    GridPath< n, Tc, Tp >& path = paths[j];
    Cell< n, Tc >& cell = path.cells.back();

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
    Edge<Cell<n,Tc>, GridPath<n,Tc> >* oldEdge = from.Find(*to, !fwd);
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
      typename vector<Cell<n, Tc> >::iterator cit = path.cells.begin();
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

    if (trackEdgesThroughCell) {
      // iterate through path and add the new edge to the list
      // of edges passing through each cell of the path
      for (typename vector< Cell< n, Tc > >::iterator cit = path.cells.begin();
           cit != path.cells.end();
           ++cit) {
        int id = grid.Id(cit->c);
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
  }

  if (fwd)
    from.succExpanded = true;
  else
    from.predExpanded = true;

  return true;
}

template < int n, class Tc, class Tp >
GridSearch< n, Tc, Tp >::~GridSearch() {
  typename std::map< int, CellEdge* >::iterator ei;
  typename std::map< int, CellVertex* >::iterator vi;

  for (ei = graph.edges.begin(); ei != graph.edges.end(); ++ei) {
    delete ei->second;
  }
  for (vi = graph.vertices.begin(); vi != graph.vertices.end(); ++vi) {
    delete vi->second;
  }

  delete[] vertexMap;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::SetStart(const Vectornd& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetStart: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);
  //    std::cout << "id=" << id << std::endl;
  assert(id >= 0 && id < grid.nc);

  // cell must exist
  Cell< n, Tc >* cell = grid.cells[id];
  if (!cell) {
    std::cout << "[W] GridSearch:SetStart: cell does not exist x="
              << x.transpose() << std::endl;
    return false;
  }

  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  //    Expand(*vertexMap[id]);

  Search< Cell< n, Tc >, GridPath< n, Tc, Tp > >::SetStart(*vertexMap[id]);

  return true;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::SetGoal(const Vectornd& x) {
  if (!grid.Valid(x)) {
    std::cout << "[W] GridSearch:SetGoal: invalid x=" << x.transpose()
              << std::endl;
    return false;
  }

  int id = grid.Id(x);
  assert(id >= 0 && id < grid.nc);

  Cell< n, Tc >* cell = grid.cells[id];
  if (!cell) {
    std::cout << "[W] GridSearch:SetGoal: cell does not exist x="
              << x.transpose() << std::endl;
    return false;
  }

  // if it's not added previously add it
  if (!vertexMap[id]) {
    vertexMap[id] = new CellVertex(*cell);
    graph.AddVertex(*vertexMap[id]);
  }

  Search< Cell< n, Tc >, GridPath< n, Tc, Tp > >::SetGoal(*vertexMap[id]);

  return true;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::RemoveCell(const Vectornd& x) {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return false;

  if (grid.cells[id]) {
    delete grid.cells[id];
    grid.cells[id] = 0;
  } else {
    return false;
  }

  CellVertex* v = vertexMap[id];
  if (v) {
    graph.RemoveVertex(*v);
  }

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

  return true;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::Plan(GridPath< n, Tc, Tp >& path,
                                   bool removeDuplicateCells) {
  path.cells.clear();
  path.cost = 0;

  std::vector< CellEdge* > edgePath;

  if (Search< Cell< n, Tc >, GridPath< n, Tc, Tp > >::Plan(edgePath) < 0)
    return false;

  typename vector< CellEdge* >::iterator it;
  for (it = edgePath.begin(); it != edgePath.end(); ++it) {
    CellEdge* edge = *it;
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
    typename vector< Cell< n, Tc > >::iterator cit;
    vector< Cell< n, Tc > > cells;
    cells.push_back(path.cells.front());
    for (cit = path.cells.begin(); (cit + 1) != path.cells.end(); ++cit) {
      if ((cit->c - (cit + 1)->c).norm() > 1e-10)
        cells.push_back(*(cit + 1));
    }
    path.cells = cells;
  }
  // now add the very last one
  //    path.cells.push_back(edgePath.back()->data.cells.back());

  return true;
}

template < int n, class Tc, class Tp >
void GridSearch< n, Tc, Tp >::OptPath(const GridPath< n, Tc, Tp >& path,
                                      GridPath< n, Tc, Tp >& optPath,
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
  typename vector< Cell< n, Tc > >::const_iterator it0;
  it0 = path.cells.begin();
  typename vector< Cell< n, Tc > >::const_iterator it1;
  it1 = it0 + 1;

  Vectornd x0 = it0->c;
  Vectornd x1 = it1->c;

  Vectornd dx0 = x1 - x0;
  double dn = dx0.norm();
  dx0 /= dn;

  optPath.cells.push_back(path.cells[0]);

  if (traceStep <= 0)
    traceStep = grid.cs.norm();

  for (; it1 != path.cells.end() - 1; ++it1) {
    x1 = it1->c;
    Vectornd x2 = (it1 + 1)->c;
    Vectornd dx1 = x2 - x0;
    dn = dx1.norm();
    assert(dn > 1e-12);
    dx1 /= dn;

    if ((dx0 - dx1).norm() > 1e-16) {
      dn = (x0 - x2).norm();
      for (double d = traceStep; d < dn; d += traceStep) {
        Vectornd x = x0 + dx1 * d;
        int id = grid.Id(x);
        assert(id >= 0 && id < grid.nc);
        if (!grid.cells[id] || grid.cells[id]->cost > freeCost) {
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

template < int n, class Tc, class Tp >
void GridSearch< n, Tc, Tp >::SplinePath(const GridPath< n, Tc, Tp >& path,
                                         std::vector< Vectornd >& splinePath,
                                         /*GridPath<n,Tc> &splineCells, */
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
    Vectornd pti;
    for (int j = 0; j < n; j++) {
      pti(j) = sps[j][i * traceStep];
    }
    splinePath[i] = pti;

    // std::cout << i*traceStep << std::endl;
    // std::cout << pti.transpose() << std::endl;
    /*
    Cell<n,Tc>* cell = grid.Get(splinePath[i]);
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
    */
  }
}

template < int n, class Tc, class Tp >
double GridSearch< n, Tc, Tp >::GetCost(const Vectornd& x) const {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return 0;

  Cell< n, Tc >* cell = grid.cells[id];
  if (!cell) {
    std::cout << "[W] GridSearch::GetCost: no cell at position "
              << x.transpose() << std::endl;
    return 0;
  }
  return cell->cost;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::SetCost(const Vectornd& x, double cost) {
  int id = grid.Id(x);
  if (id < 0 || id >= grid.nc)
    return false;

  Cell< n, Tc >* cell = grid.cells[id];
  if (!cell) {
    // std::cout << "[W] GridSearch::SetCost: no cell at position " <<
    // x.transpose() << std::endl;
    return false;
  }

  //    std::cout << "[W] GridSearch::SetCost: checking position " <<
  //    x.transpose() << std::endl;

  // if the cost has not changed simply return
  if (Search< Cell< n, Tc >, GridPath< n, Tc, Tp > >::Eq(cost, cell->cost))
    return false;

  cell->cost = cost;

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

  vertex->data.cost = cost;

  // fix edges are not tracked, then at best we can modify the once incoming and
  // outgoing from the cell
  // if edges are tracked this is already done above
  if (!trackEdgesThroughCell) {
    for (typename std::map< int, CellEdge* >::iterator ein = vertex->in.begin();
         ein != vertex->in.end();
         ein++) {
      this->ChangeCost(*ein->second,
                       this->cost.Real(ein->second->from->data, vertex->data));
    }

    for (typename std::map< int, CellEdge* >::iterator eout =
             vertex->out.begin();
         eout != vertex->out.end();
         eout++) {
      this->ChangeCost(*eout->second,
                       this->cost.Real(vertex->data, eout->second->to->data));
    }
  }
  return true;
}

template < int n, class Tc, class Tp >
bool GridSearch< n, Tc, Tp >::AddEdge(const Vectornd& x1, const Vectornd& x2) {
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
