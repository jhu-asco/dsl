// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "gridconnectivity.h"
#include "grid.h"
#include "spline.h"

namespace dsl {

/**
 * The most basic grid connectivity: each cell is connected with "lines" to
 * other cells, and each line has an associate cost.
 *
 * DataT defines the cost of going through a point and
 * should be a primitives type bool,int,float,double, etc...
 *
 * Author: Marin Kobilarov
 */
template < class PointT, class DataT = EmptyData >
class LineConnectivity : public GridConnectivity< PointT, DataT, PointT > {
public:
  using CellT = Cell< PointT, DataT >;
  using PathT = GridPath< PointT, DataT, PointT >;

  using PrimitiveT = std::tuple< Cell< PointT, DataT >*, PointT, double >;

  /**
   * Initialize connectivity using a grid
   * @param grid the grid
   */
  LineConnectivity(const Grid< PointT, DataT >& grid);

  /**
   * Initialize connectivity using a grid, lines, and costs
   * @param grid the grid
   * @param lines the lines
   * @param costs the costs
   */
  LineConnectivity(const Grid< PointT, DataT >& grid,
                   const std::vector< PointT >& lines,
                   const std::vector< double >& costs);

  // the data stored in each cell is just the occupancy probability
  // if the occupancy probability is less than 0.5 we return free
  bool free(const DataT& cost) const override {
    return cost < 0.5;
  }

  bool operator()(const Cell< PointT, DataT >& from,
                  std::vector< PrimitiveT >& paths,
                  bool fwd = true) const override;

  /**
   * Experimental path "straightening" function
   * @param path original path
   * @param opt_path optimized path
   * @param freeCost (anything above freeCost is considered an obstacle through
   * the path cannoth pass)
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is -1 which means the
   * internally the cell size is used)
   */
  void straightPath(const PathT& path,
                    PathT& opt_path,
                    double freeCost = 1e-3,
                    double traceStep = -1.0) const;

  /**
   * Experimental path "smoothing" function
   * @param path original path
   * @param splinePath spline path
   * @param traceStep step with which to trace the path during optimization
   * (should be comparable to the cell size, by defaut is 0.1)
   */
  void splinePath(const PathT& path,
                  std::vector< PointT >& splinePath,
                  // GridPath<n,CellData> &splineCells,
                  double traceStep = 0.1) const;

  const Grid< PointT, DataT >& grid; ///< grid
  std::vector< PointT >
      lines; ///< line vectors (directions) connecting to other cells
  std::vector< double > costs; ///< cost along each direction
};

template < class PointT, class DataT >
LineConnectivity< PointT, DataT >::LineConnectivity(
    const Grid< PointT, DataT >& grid)
  : grid(grid) {}

template < class PointT, class DataT >
LineConnectivity< PointT, DataT >::LineConnectivity(
    const Grid< PointT, DataT >& grid,
    const std::vector< PointT >& lines,
    const std::vector< double >& costs)
  : grid(grid), lines(lines), costs(costs) {}

template < class PointT, class DataT >
bool LineConnectivity< PointT, DataT >::
    operator()(const Cell< PointT, DataT >& from,
               std::vector< PrimitiveT >& paths,
               bool fwd) const {
  paths.clear();
  for (size_t i = 0; i < lines.size(); ++i) {
    PointT x = from.center;
    if (fwd)
      x += lines[i];
    else
      x -= lines[i];

    if (!grid.valid(x))
      continue;

    CellT* cell = grid.data(x, false);
    if (!cell)
      continue;

    if (!free(cell->data)) // if obstacle
      continue;

    paths.push_back(std::make_tuple(cell, lines[i], costs[i]));
  }
  return true;
}

template < class PointT, class DataT >
void LineConnectivity< PointT, DataT >::straightPath(const PathT& path,
                                                     PathT& opt_path,
                                                     double freeCost,
                                                     double traceStep) const {
  double len = 0;

  opt_path.cells.clear();
  opt_path.cost = 0;

  if (path.cells.size() == 2) {
    opt_path.cells = path.cells;
    opt_path.cost = path.cost;
    return;
  }
  auto it0 = path.cells.begin();
  auto it1 = it0 + 1;

  PointT x0 = it0->center;
  PointT x1 = it1->center;

  PointT dx0 = x1 - x0;
  double dn = dx0.norm();
  dx0 /= dn;

  opt_path.cells.push_back(path.cells[0]);

  if (traceStep <= 0)
    traceStep = grid.cs.norm();

  for (; it1 != path.cells.end() - 1; ++it1) {
    x1 = it1->center;
    PointT x2 = (it1 + 1)->center;
    PointT dx1 = x2 - x0;
    dn = dx1.norm();
    assert(dn > 1e-12);
    dx1 /= dn;

    if ((dx0 - dx1).norm() > 1e-16) {
      dn = (x0 - x2).norm();
      for (double d = traceStep; d < dn; d += traceStep) {
        PointT x = x0 + dx1 * d;
        int id = grid.computeId(x);
        assert(id >= 0 && id < grid.size);
        if (!grid.values[id] || !free(grid.values[id]->data)) {
          opt_path.cells.push_back(*it1);
          x0 = x1;
          break;
        }
      }
      dx0 = x2 - x0;
      dn = dx0.norm();
      dx0 /= dn;
    }
  }

  opt_path.cells.push_back(path.cells.back());
}

template < class PointT, class DataT >
void LineConnectivity< PointT, DataT >::splinePath(
    const GridPath< PointT, DataT, PointT >& path,
    std::vector< PointT >& splinePath,
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

  for (int i = 0; i < count; i++) {
    PointT pti;
    for (int j = 0; j < n; j++) {
      pti(j) = sps[j][i * traceStep];
    }
    splinePath[i] = pti;

    //  TODO(comment this out)
    Cell< PointT, DataT >* cell = grid.data(splinePath[i]);
    if(!cell)
    {
      std::cout << "dsl::splinePath: Cell does not exist" << std::endl;
      return;
    }
    /*
    splineCells.cells[i] = *(grid.data(splinePath[i]));
    if(i > 0)
    {
      splineCells.cost +=
    (splineCells.cells[i].c-splineCells.cells[i-1].c).norm();
    }
    */
  }
}
}
