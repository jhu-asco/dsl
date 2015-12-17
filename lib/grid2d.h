// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID2D_H
#define DSL_GRID2D_H

#include "grid.h"

namespace dsl {
  
  class Grid2d : public Grid<2> {
  public:
    
    Grid2d(int width, int height, double *map, 
           double sx, double sy, double costScale,
           double maxCost = 1);
  };
}

#endif
