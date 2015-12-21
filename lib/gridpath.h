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
   *  Path containing a list of grid points
   */
  template<int n>
    class GridPath {
  public:
    
  GridPath() : len(0) {
    };

    std::vector<Cell<n> > cells;  ///< list of cells along path
    double len;                   ///< length of path (sum of edge b/n cells)

  };
}

#endif
