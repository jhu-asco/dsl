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
    class LineConnectivity : public GridConnectivity< PointType, DataType > {
public:

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
                   const std::vector< Cell<PointType, DataType> >& lines,
                   const std::vector< double >& costs);

  virtual bool Free(const DataType &cost) const override { return cost < 0.5; }

  virtual bool operator()(const Cell<PointType, DataType>& from,
                          std::vector< GridPath< PointType, DataType> >& paths,
                          bool fwd = true) const;

  const Grid<  PointType, DataType >& grid; ///< grid

  std::vector< PointType > lines; ///< line vectors (directions) connecting to other cells
  std::vector< double > costs; ///< cost along each direction
};

template < class PointType, class DataType >
    LineConnectivity< PointType, DataType >::LineConnectivity(const Grid< PointType, DataType >& grid)
  : grid(grid) {}

template < class PointType, class DataType >
LineConnectivity< PointType, DataType >::LineConnectivity(const Grid<  PointType, DataType >& grid,
                                                          const std::vector< Cell< PointType, DataType> >& lines,
                                                          const std::vector< double >& costs)
  : grid(grid), lines(lines), costs(costs) {}

template < class PointType, class DataType >
    bool LineConnectivity< PointType, DataType >::operator()(const Cell< PointType, DataType >& from,
                                                             std::vector< GridPath< PointType, DataType > >& paths,
                                                             bool fwd) const {
  paths.clear();
  for (size_t i = 0; i < lines.size(); ++i) {
    GridPath< PointType, DataType > path;
    path.cells.push_back(from);

    PointType x = from.c;
    if (fwd)
      x += lines[i];
    else
      x -= lines[i];

    if (!grid.Valid(x))
      continue;

    auto cell = grid.Get(x, false);
    if (!cell)
      continue;

    if (!Free(cell->data)) // if obstacle
      continue;
    
    //    path.cells.push_back(Cell<PointType, DataType>(0, x, *data));
    path.cells.push_back(*cell);
    //    path.cost = (1 + (from.data + data) / 2) * this->costs[i];
    path.cost = this->costs[i];
    path.fwd = fwd;
    paths.push_back(path);
  }
  return true;
}
}

#endif
