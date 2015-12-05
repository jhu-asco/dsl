// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDCOST3D_H
#define DSL_GRIDCOST3D_H

#include "cost.h"
#include "cell3d.h"

namespace dsl {
  /**
   * 3D Grid cost interface.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  class GridCost3D : public Cost<Cell3d> {
  public:
    double Heur(const Cell3d &va, const Cell3d &vb) const;       
    double Real(const Cell3d &va, const Cell3d &vb) const;    
  };
}

#endif
