// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "gridcost.h"

namespace dsl {

/**
 * Traversability cost, defining cost between two cells,
 * the cell data is interpreted as traversability (0=free)
 *
 * Author: Marin Kobilarov -- Copyright (C) 2004
 */
  template < class PointType, class DataType >
      class TraversabilityCost : public GridCost< PointType, DataType > {
public:

  using CellType = Cell<PointType, DataType>;
    
  double real(const CellType& a, const CellType& b) const {
    // default real cost is euclidean distance + average cell cost multiplied by
    // Euclidean distance
    //return (a.c - b.c).norm();
    return (1 + (a.data + b.data) / 2) * (a.c - b.c).norm();
  }

  double heur(const CellType& a, const CellType& b) const {
    /*
    double dx = std::abs(a.c[0] - b.c[0]);
    double dy = std::abs(a.c[1] - b.c[1]);
    if (dx > dy)
      return dx;
    else
      return dy;
    */
      
    // default Heuristic cost is the Euclidean distance
    //return (1 + (a.data + b.data) / 2) * (a.c - b.c).norm();
    return  .99*(a.c - b.c).norm();
  }
};
}
