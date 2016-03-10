// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID_H
#define DSL_GRID_H

#include "cell.h"
#include <Eigen/Dense>

namespace dsl {

struct EmptyData {};

/**
 * An n-dimenensional grid consisting of abstract "cells", or elements
 * identified by a set of coordinates of type PointType, each cell
 * containing data of type DataType.
 * A grid provides instant access to the elements by maintaining
 * an n-dimensional array of pointers to cells. Cells that are empty,
 * e.g. that are inside obstacles, correspond to null pointers.
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 */
template < class PointType, class DataType = EmptyData>
 struct Grid {

 using Vectorni =  Eigen::Matrix< int, PointType::SizeAtCompileTime, 1 >;
 using TypedCell = Cell<PointType, DataType>;
   
  /**
   * Initialize the grid using state lower bound, state upper bound, the number
   * of grid cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of grid cells per each dimension
   */
  Grid(const PointType& xlb, const PointType& xub, const Vectorni& gs)
      : n(xlb.size()), xlb(xlb), xub(xub), gs(gs) {
    ds = xub - xlb;
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      nc *= gs[i]; // total number of cells
      cs[i] = (xub[i] - xlb[i]) / gs[i];
    }    
    cells = new TypedCell*[nc];
    memset(cells, 0, nc * sizeof(TypedCell*)); // initialize all of them nil
  }

    /**
   * Initialize the map using state lower bound, state upper bound, the number
   * of map cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param cs cell dimensions
   */
  Grid(const PointType& xlb, const PointType& xub, const PointType& cs)
      : n(xlb.size()), xlb(xlb), xub(xub), cs(cs) {
    ds = xub - xlb;
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] < xub[i]);
      assert(cs[i] > 0);
      gs[i] = floor((xub[i] - xlb[i]) / cs[i]);
      nc *= gs[i]; // total number of cells
    }

    cells = new TypedCell*[nc];
    memset(cells, 0, nc * sizeof(TypedCell*)); // initialize all of them nil
  }

  Grid(const Grid &grid) : n(grid.n), xlb(grid.xlb), xub(grid.xub), ds(grid.ds), cs(grid.cs), gs(grid.gs), nc(grid.nc) {
     cells = new TypedCell*[nc];
     memcpy(cells, grid.cells, nc * sizeof(TypedCell*)); // initialize all of them nil
   }

  virtual ~Grid() {
    delete[] cells;
  }

  /**
   * Check if point x is within grid bounds
   * @param x point
   * @return true if within bounds
   */
  virtual bool Valid(const PointType& x, double eps = 1e-10) const {
    for (int i = 0; i < x.size(); ++i) {
      if (x[i] < xlb[i] + eps)
        return false;
      if (x[i] > xub[i] - eps)
        return false;
    }
    return true;
  }

  /**
   * Get an id of point x useful for direct lookup in the grid array
   * @param x point
   * @return a computed id
   */
  int Id(const PointType& x) const {

       // unroll loop for n=1,2,3,4 for efficiency
    if (n==1)
      return Index(x, 0);

    if (n==2)
      return Index(x, 0) + gs[0]*Index(x, 1);

    if (n==3)
      return Index(x, 0) + gs[0]*Index(x, 1) + gs[0]*gs[1]*Index(x, 2);

    if (n==4) {
      int cum = gs[0]*gs[1];
      return Index(x, 0) + gs[0]*Index(x, 1) + cum*Index(x, 2) + cum*gs[2]*Index(x, 3);
    }

    int cum = 1; // cumulative offset for next dimension
    int id = 0;
    for (int i = 0; i < x.size(); ++i) {
      // index of i-th dimension
      int ind = floor((x[i] - xlb[i]) / (xub[i] - xlb[i]) * gs[i]);
      id += cum * ind;
      if (i < n - 1)
        cum *= gs[i];
    }
    return id;
  }

  /**
   * Get the grid index of i-th dimension of point x
   * @param x point
   * @param i coordinate index
   * @return index into cell array
   */
  int Index(const PointType& x, int i) const {
    return floor((x[i] - xlb[i]) / (xub[i] - xlb[i]) * gs[i]);
  }

   
  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return pointer to a cell or 0 if cell does not exist
   */
   TypedCell* Get(const PointType& x, bool checkValid = true) const {
    if (checkValid)
      if (!Valid(x))
        return 0;
    return Get(Id(x));
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
   TypedCell* Get(int id) const {
    assert(id >= 0);
    if (id >= nc)
      return 0;
    return cells[id];
  }

  int n;      ///< grid dimension
  
  PointType xlb; ///< state lower bound
  PointType xub; ///< state upper bound
  PointType ds;  ///< dimensions (ds=xub-xlb)
  PointType cs;  ///< cell length size per dimension
  Vectorni gs;  ///< number of cells per dimension

  int nc = 0;        ///< number of cells in grid
  TypedCell** cells = nullptr; ///< array of cell data
};
}

#endif
