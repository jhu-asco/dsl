// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDPATH_H
#define DSL_GRIDPATH_H

#include "cell.h"
#include <vector>

namespace dsl {

/**
 *  A generic grid "path" represented by:
 * 1) a list of cells it passes through
 * 2) a generic data type that is application-specific
 * 3) the total cost (e.g. distance,time,etc...) of the path
 *
 * The cells containt generic data type CellData, while the grid path contains generic
 * data type PathData, which by default is a vector of CellData's.  PathData is useful for
 *representing
 * e.g. a continuous trajectory corresponding to the discrete sequence of cells.
 *
 */
 template <class PointType, class DataType>
class GridPath {
public:
  GridPath() : cost(0), fwd(0) {};

   std::vector< Cell<PointType, DataType> > cells; ///< list of cells along path

   // PathData data; ///< generic data stored along path; typically this represents a
  /// path/trajectory passing through the cells

  double cost; ///< cost of path

  bool fwd; ///< generated forward or backward
};
}

#endif
