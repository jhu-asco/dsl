// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARGRID_H
#define DSL_CARGRID_H

#include "grid.h"

namespace dsl {
  
  class CarGrid : public Grid<3> {
  public:
    
    CarGrid(int width, int height, double *map, 
            double sx, double sy, double sa, double costScale,
            double maxCost = 1);

    double maxCost; ///< any cell cost above maxCost is considered obstacle and not added to the graph
  };
}

#endif
