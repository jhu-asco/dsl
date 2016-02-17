// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LINECONNECTIVITY_H
#define DSL_LINECONNECTIVITY_H

#include "gridconnectivity.h"
#include "grid.h"
#include "spline.h"

namespace dsl {

/**
 * The most basic grid connectivity: each cell is connected with "lines" to
 * other cells, and each line has an associate cost.
 *
 * DataType defines the cost of going through a point and
 * should be a primitives type bool,int,float,double, etc...
 *
 * Author: Marin Kobilarov
 */
template < class PointType, class DataType = EmptyData>
    class LineConnectivity : public GridConnectivity< PointType, DataType, PointType> {
public:

using TypedCell = Cell<PointType, DataType>;

  /**
   * Initialize connectivity using a grid
   * @param grid the grid
   */
  LineConnectivity(const Grid< PointType, DataType >& grid);

  /**
   * Initialize connectivity using a grid, lines, and costs
   * @param grid the grid
   * @param lines the lines
   * @param costs the costs
   */
  LineConnectivity(const Grid<PointType, DataType>& grid,
                   const std::vector<PointType>& lines,
                   const std::vector< double >& costs);

  virtual bool Free(const DataType &cost) const override { return cost < 0.5; }

  virtual bool operator()(const Cell<PointType, DataType>& from,
                          std::vector< std::tuple<Cell<PointType, DataType>*, PointType, double> >& paths,
                          bool fwd = true) const;

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
  void OptPath(const GridPath< PointType, DataType, PointType >& path,
               GridPath< PointType, DataType, PointType >& optPath,
               double freeCost = 1e-3,
               double traceStep = -1.0) const;

   

  /**
   * Experimental path "smoothing" function
   * @param path original path
   * @param splinePath spline path
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is 0.1)
   */
  void SplinePath(const GridPath< PointType, DataType, PointType >& path,
                  std::vector< PointType >& splinePath,
                  // GridPath<n,CellData> &splineCells,
                  double traceStep = 0.1) const;


  const Grid<  PointType, DataType >& grid; ///< grid
  std::vector< PointType > lines; ///< line vectors (directions) connecting to other cells
  std::vector< double > costs; ///< cost along each direction
};

template < class PointType, class DataType >
    LineConnectivity< PointType, DataType >::LineConnectivity(const Grid< PointType, DataType >& grid)
  : grid(grid) {}

template < class PointType, class DataType >
LineConnectivity< PointType, DataType >::LineConnectivity(const Grid<  PointType, DataType >& grid,
                                                          const std::vector< PointType >& lines,
                                                          const std::vector< double >& costs)
  : grid(grid), lines(lines), costs(costs) {}

template < class PointType, class DataType >
    bool LineConnectivity< PointType, DataType >::operator()(const Cell< PointType, DataType >& from,
                                                             std::vector< std::tuple<Cell<PointType, DataType>*, PointType, double> >& paths,
                                                             bool fwd) const {
  paths.clear();
  for (size_t i = 0; i < lines.size(); ++i) {
    // GridPath< PointType, DataType > path;
    //  path.cells.push_back(from);

    PointType x = from.c;
    if (fwd)
      x += lines[i];
    else
      x -= lines[i];

    if (!grid.Valid(x))
      continue;

    TypedCell *cell = grid.Get(x, false);
    if (!cell)
      continue;

    if (!Free(cell->data)) // if obstacle
      continue;

    paths.push_back(std::make_tuple(cell, lines[i], costs[i]));
    
    //    path.cells.push_back(*cell);
    //    path.cost = this->costs[i];
    //    path.fwd = fwd;
    //    paths.push_back(path);
  }
  return true;
}



template < class PointType, class DataType >
    void LineConnectivity< PointType, DataType>::OptPath(
        const GridPath< PointType, DataType, PointType >& path,
        GridPath< PointType, DataType, PointType >& optPath,
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
        if (!grid.cells[id] || !Free(grid.cells[id]->data)) {
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


template < class PointType, class DataType >
    void LineConnectivity< PointType, DataType >::SplinePath(const GridPath< PointType, DataType, PointType>& path,
                                                             std::vector< PointType >& splinePath,
                                                             double traceStep) const {
  std::vector< double > steps(path.cells.size());
  int n = lines[0].size();
  std::vector< std::vector< double > > pts(n, std::vector< double >(path.cells.size()));
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
    Cell<PointType, DataType>* cell = grid.Get(splinePath[i]);
    if(!cell)
    {
      std::cout << "dsl::SplinePath: Cell does not exist" << std::endl;
      return;
    }
    /*
    splineCells.cells[i] = *(grid.Get(splinePath[i]));
    if(i > 0)
    {
      splineCells.cost +=
    (splineCells.cells[i].c-splineCells.cells[i-1].c).norm();
    }
    */
  }
}



}

#endif
