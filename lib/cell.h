// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CELL_H
#define DSL_CELL_H

#include <Eigen/Dense>

namespace dsl {

using namespace Eigen;

/**
 * A basic grid cell defining a rectangular region in n-dimension space
 * The cell can also contain a generic data of type T, which by default is just
 * a point in the cell.
 */
template < int n, class Tc = Matrix< double, n, 1 > >
class Cell {
  typedef Matrix< double, n, 1 > Vectornd;

public:
  Cell() : cost(0) {}

  /**
   * Initialize a cell using its center and cost
   * @param c center
   * @param cost cost
   */
  Cell(const Vectornd& c, double cost = 0) : c(c), cost(cost) {}

  /**
   * Initialize a cell using its center and cost
   * @param c center
   * @param r half-distance of cell size (i.e. max radius of sphere fitting
   * inside cell)
   * @param cost cost
   */
  Cell(const Vectornd& c, const Vectornd& r, double cost = 0)
    : c(c), r(r), cost(cost) {}

  Cell(const Cell< n, Tc >& cell)
    : c(cell.c), r(cell.r), cost(cell.cost), data(cell.data) {}

  /**
   * Is a point x inside the cell
   * @param x point
   * @return true if inside
   */
  bool Inside(const Vectornd& x) {
    for (int i = 0; i < x.size(); ++i)
      if (x[i] < c[i] - r[i] || x[i] > c[i] + r[i])
        return false;
    return true;
  }

  Vectornd c; ///< center of cell
  Vectornd r; ///< half-distance of each cell side

  double cost; ///< cost of cell (typically this is 0 if unoccupied)

  Tc data; ///< generic data stored in the cell (e.g. could store the actual
  /// position of a system which might not coincide with the center of
  /// the cell, but more generally can store extra algorithm-specific
  /// data)
};
}

#endif
