// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "cost.h"
#include "gridpath.h"

namespace dsl {

/**
 * Grid cost interface, defining cost between two cells
 *
 * Author: Marin Kobilarov -- Copyright (C) 2004
 */
  template < class PointType, class DataType >
class GridCost : public Cost< Cell<PointType, DataType> > {
public:

  using GridVertexData = Cell<PointType, DataType>;

  double real(const GridVertexData& a, const GridVertexData& b) const {
    // default real cost is euclidean distance
    return (a.centr - b.centr).norm();
    //  another version adds a term for the average cell cost multiplied by
    // Euclidean distance, as specified below:
    // return (1 + (a.cost + b.cost) / 2) * (a.centr - b.centr).norm();
  }

  double heur(const GridVertexData& a, const GridVertexData& b) const {
    // default Heuristic cost is the Euclidean distance
    return (a.centr - b.centr).norm();
  }
};
}
