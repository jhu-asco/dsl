// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GRIDCONNECTIVITY_H_
#define DSL_LIB_GRIDCONNECTIVITY_H_

#include "gridpath.h"
#include <tuple>
#include <memory>

namespace dsl {

  template < class PointType, class DataType, class ConnectionType>
class GridConnectivity {
public:
    using TypedCell = Cell<PointType, DataType>;
    using TypedCellPtr = std::shared_ptr<Cell<PointType, DataType> >;

  /**
   * Connectivity operator providing primitive paths from a given vertex. This
   * function should generate connections from the given vertex and store them
   * in
   * the paths object.
   * @param from starting vertex
   * @param paths generated paths
   * @param fwd true if generated forward in time
   * @return true on success
   */
  virtual bool operator()(const TypedCell& from,
                          std::vector< std::tuple<TypedCellPtr, ConnectionType, double> >& toCells,
                          bool fwd = true) const = 0;

  virtual bool Free(const DataType &data) const = 0;

};
}

#endif
