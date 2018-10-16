// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_TRAVERSABILITYCOST_H
#define DSL_TRAVERSABILITYCOST_H

#include "gridcost.h"

namespace dsl {

/**
 * Traversability cost, defining cost between two cells,
 * the cell data is interpreted as traversability (0=free)
 *
 * Author: Marin Kobilarov -- Copyright (C) 2004
 */
template < class PointT, class DataT >
class TraversabilityCost : public GridCost< PointT, DataT > {
public:
  using CellT = Cell< PointT, DataT >;

  double real(const CellT& a, const CellT& b) const {
    // default real cost is euclidean distance + average cell cost multiplied by
    // Euclidean distance
    // return (a.center - b.center).norm();
    return (1 + (a.data + b.data) / 2) * (a.center - b.center).norm();
  }

  double heur(const CellT& a, const CellT& b) const {
    double dx = std::abs(a.center[0] - b.center[0]);
    double dy = std::abs(a.center[1] - b.center[1]);
    if (dx > dy)
      return dx;
    else
      return dy;

    // default Heuristic cost is the Euclidean distance
    // return (1 + (a.data + b.data) / 2) * (a.center - b.center).norm();
    return .99 * (a.center - b.center).norm();
  }
};
}

#endif
