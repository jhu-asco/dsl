// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID3DCONNECTIVITY_H
#define DSL_GRID3DCONNECTIVITY_H

#include "lineconnectivity.h"

namespace dsl {
  
  
  class Grid3dConnectivity : public LineConnectivity<3> {
  public:
    Grid3dConnectivity(const Grid<3> &grid);
  };
}

#endif
