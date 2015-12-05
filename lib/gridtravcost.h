// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDTRAVCOST_H
#define DSL_GRIDTRAVCOST_H

#include "gridcost.h"

/**
 *  Calculates the edge cost as the height gradient between two vertices 
 *
 *  Author: Matt Sheckells (c) 2015 msheckells(at)jhu.edu
 */


namespace dsl {

  class GridTravCost : public GridCost
  {
  public:
    
    double Real(const Cell2d &va, const Cell2d &vb) const;       
  };
}


#endif
