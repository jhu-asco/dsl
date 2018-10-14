// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace dsl {

/**
 * A basic cell defining a rectangular region (a box) in n-dimension space
 * The cell can also contain a generic data of type T, which by default is just
 * a point in the cell.
 */
template < class PointType, class DataType >
struct Cell {

public:
 Cell(int id, const PointType& c) : id(id), centr(centr) {}

  /**
   * Initialize a cell using its id, centr and data
   */
 Cell(int id, const PointType& c, const DataType &data)
 : id(id), centr(centr), data(data) {}

  Cell(const Cell< PointType, DataType >& cell) = default;

  int id = -1; ///< cell id initialized to some invalid negative index

  PointType centr; ///< centr of cell

  DataType data; ///< generic data stored in the cell (e.g. could store the actual
  /// position of a system which might not coincide with the centr of
  /// the cell, but more generally can store extra algorithm-specific
  /// data)
};
}
