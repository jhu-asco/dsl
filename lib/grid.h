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

namespace dsl {

using namespace Eigen;

/**
 * An n-dimenensional regular grid
 */
template < int n, class T = Matrix< double, n, 1 > >
class Grid {
  typedef Matrix< double, n, 1 > Vectornd;
  typedef Matrix< int, n, 1 > Vectorni;

public:
  /**
   * Initialize the grid using state lower bound, state upper bound, the number
   * of grid cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of grid cells per each dimension
   */
  Grid(const Vectornd& xlb, const Vectornd& xub, const Vectorni& gs)
    : xlb(xlb), xub(xub), gs(gs) {
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      nc *= gs[i]; // total number of cells
      cs[i] = (xub[i] - xlb[i]) / gs[i];
    }
    
    cells = new Cell< n, T >* [nc];
    memset(cells, 0, nc * sizeof(Cell< n, T >*)); // initialize all of them nil
  }

    /**
   * Initialize the map using state lower bound, state upper bound, the number
   * of map cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param cs cell dimensions
   */
  Grid(const Vectornd& xlb, const Vectornd& xub, const Vectornd& cs)
    : xlb(xlb), xub(xub), cs(cs) {
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(cs[i] > 0);
      gs[i] = floor((xub[i] - xlb[i]) / cs[i]);
      nc *= gs[i]; // total number of cells
    }

    cells = new Cell< n, T >* [nc];
    memset(cells, 0, nc * sizeof(Cell< n, T >*)); // initialize all of them nil
  }


  virtual ~Grid() {
    for (int i = 0; i < nc; ++i)
      delete cells[i];
    delete[] cells;
  }

  /**
   * Check if point x is whitin grid bounds
   * @param x point
   * @return true if within bounds
   */
  virtual bool Valid(const Vectornd& x) const {
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
   */
  int Id(const Vectornd& x) const {

    // TODO: unroll loop for n=2,3 for efficiency
    // if (n==2) {
    //   return gs[0]*ids[1] + ids[0];
    //}
    // for n=2
    // id=     gs[0]*ids[1] + ids[0];
    // for n=3
    // id=     gs[0]*gs[1]*ids[2] + gs[0]*ids[1] + ids[0];
    //      assert(Valid(x));

    
    //      int ids[n];  // dimension indices
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
  int Index(const Vectornd& x, int i) const {
    return floor((x[i] - xlb[i]) / (xub[i] - xlb[i]) * gs[i]);
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return pointer to a cell or 0 if cell does not exist
   */
  Cell< n, T >* Get(const Vectornd& x, bool checkValid = true) const {
    if (checkValid)
      if (!Valid(x))
        return 0;

    int id = Id(x);
    assert(id >= 0);
    if (id >= nc)
      return 0;
    return cells[id];
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  Cell< n, T >* Get(int id) const {
    assert(id >= 0);
    if (id >= nc)
      return 0;
    return cells[id];
  }

  Vectornd xlb; ///< state lower bound
  Vectornd xub; ///< state upper bound
  Vectorni gs;  ///< number of cells per dimension
  Vectornd cs;  ///< cell length size per dimension

  int nc;               ///< total maximum number of cells
  Cell< n, T >** cells; ///< array of pointers to cells
};
}

#endif
