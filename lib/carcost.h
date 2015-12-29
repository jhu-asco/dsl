// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCOST_H
#define DSL_CARCOST_H

#include "gridcost.h"

namespace dsl {

  typedef Cell<3, Matrix3d> SE2Cell;

  /**
   * Grid cost interface.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  class CarCost : public GridCost<3, Matrix3d> {
  public:
    double Heur(const SE2Cell &a, const SE2Cell &b) const;       
    double Real(const SE2Cell &a, const SE2Cell &b) const;    
  };
}

#endif
