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

  /**
   * Grid cost interface.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  class CarCost : public GridCost<3 > {
  public:
    double Heur(const Cell<3> &a, const Cell<3> &b) const;       
    double Real(const Cell<3> &a, const Cell<3> &b) const;    
  };
}

#endif
