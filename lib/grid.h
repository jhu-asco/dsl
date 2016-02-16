// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID_H
#define DSL_GRID_H

//#include "lattice.h"
#include "cell.h"
#include <Eigen/Dense>

namespace dsl {

struct EmptyData {};

/**
 * An n-dimenensional regular grid
 */
template < class PointType = Eigen::VectorXd, class DataType = EmptyData>
 struct Grid {

 using Vectorni =  Eigen::Matrix< int, PointType::SizeAtCompileTime, 1 >;
 using TCell = Cell<PointType, DataType>;
   
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
    cells = new TCell*[nc];
    memset(cells, 0, nc * sizeof(TCell*)); // initialize all of them nil
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

    cells = new TCell*[nc];
    memset(cells, 0, nc * sizeof(TCell*)); // initialize all of them nil
  }

Grid(const Grid &grid) : n(grid.n), xlb(grid.xlb), xub(grid.xub), ds(grid.ds), cs(grid.cs), gs(grid.gs), nc(grid.nc) {
     cells = new TCell*[nc];
     memcpy(cells, grid.cells, nc * sizeof(TCell*)); // initialize all of them nil
   }

  virtual ~Grid() {
    delete[] cells;
  }

  /**
   * Check if point x is whitin grid bounds
   * @param x point
   * @return true if within bounds
   */
  virtual bool Valid(const PointType& x) const {
    for (int i = 0; i < x.size(); ++i) {
      if (x[i] < xlb[i])
        return false;
      if (x[i] > xub[i])
        return false;
    }
    return true;
  }

  /**
   * Get an id of point x useful for direct lookup in the grid array
   * @param x point
   * @return a computed id
  PointType Point(int id) const {
       // unroll loop for n=1,2,3,4 for efficiency
    return floor((x[i] - xlb[i]) / ds[i] * gs[i]);

    if (n==1)
      return PointType(xlb[0] + id*ds[i]/gs[0]);

    if (n==2)
      return (x[0] - xlb[0]) / ds[0] * gs[0]) + gs[0]*(x[1] - xlb[1]) / ds[1] * gs[1]);

    if (n==3)
      return Index(x, 0) + gs[0]*Index(x, 1) + gs[0]*gs[1]*Index(x, 2);

    if (n==4) {
      int cum = gs[0]*gs[1];
      return Index(x, 0) + gs[0]*Index(x, 1) + cum*Index(x, 2) + cum*gs[2]*Index(x, 3);
    }
  }
   */

  
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
   const TCell* Get(const PointType& x, bool checkValid = true) const {
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
   const TCell* Get(int id) const {
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

  int nc;        ///< number of cells in grid
  TCell** cells; ///< array of cell data
    
   //  DataType empty;
};
}

#endif
