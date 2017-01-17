// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_CELL_H_
#define DSL_LIB_CELL_H_

#include <memory>

namespace dsl {

/**
 * A basic cell defining a rectangular region (a box) in n-dimension space
 * The cell can also contain a generic data of type T, which by default is just
 * a point in the cell.
 */
template < class PT, class DT >
struct Cell {

using Ptr = std::shared_ptr< Cell<PT,DT> >;
using PointType = PT;
using DataType = DT;
 Cell(int id, const PointType& c) : id(id), c(c) {}
  
  /**
   * Initialize a cell using its center and cost
   * @param c center
   * @param data data
   */
 Cell(int id, const PointType& c, const DataType &data)
 : id(id), c(c), data(data) {}

  Cell(const Cell< PointType, DataType >& cell)
  : id(cell.id), c(cell.c), data(cell.data) {}

  /**
   * Is a point x inside the cell
   * @param x point
   * @return true if inside
  bool Inside(const PointType& x) {
    for (int i = 0; i < x.size(); ++i)
      if (x[i] < c[i] - r[i] || x[i] > c[i] + r[i])
        return false;
    return true;
  }
   */

  int id = -1;
  
  PointType c; ///< center of cell
  //  PointType r; ///< half-distance of each cell side

  //  double cost; ///< cost of cell (typically this is 0 if unoccupied)

  DataType data; ///< generic data stored in the cell (e.g. could store the actual
  /// position of a system which might not coincide with the center of
  /// the cell, but more generally can store extra algorithm-specific
  /// data)
};
}

#endif
