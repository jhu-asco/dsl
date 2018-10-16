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
template < class PointT, class DataT >
class GridCost : public Cost< Cell< PointT, DataT > > {
public:
  using VertexDataT = Cell< PointT, DataT >;

  double real(const VertexDataT& a, const VertexDataT& b) const {
    // default real cost is euclidean distance
    return (a.center - b.center).norm();
    //  another version adds a term for the average cell cost multiplied by
    // Euclidean distance, as specified below:
    // return (1 + (a.cost + b.cost) / 2) * (a.center - b.center).norm();
  }

  double heur(const VertexDataT& a, const VertexDataT& b) const {
    // default Heuristic cost is the Euclidean distance
    return (a.center - b.center).norm();
  }
};
}
