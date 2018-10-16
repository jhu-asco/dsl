// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "gridpath.h"
#include <tuple>

namespace dsl {

template < class PointT, class DataT, class ConnectionT >
class GridConnectivity {
public:
  using CellT = Cell< PointT, DataT >;

  /**
   * Connectivity operator providing primitive paths from a given vertex. This
   * function should generate connections from the given vertex and store them
   * in the paths object.
   * @param from starting vertex
   * @param paths generated paths
   * @param fwd true if generated forward in time
   * @return true on success
   */
  virtual bool operator()(
      const CellT& from,
      std::vector< std::tuple< CellT*, ConnectionT, double > >& to_cells,
      bool fwd = true) const = 0;

  virtual bool free(const DataT& data) const = 0;
};
}
