// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "gridtravcost.h"

using namespace std;
using namespace dsl;


double GridTravCost::Real(const Cell2d &va, const Cell2d &vb) const
{
  int dx = vb.p[0]-va.p[0];
  int dy = vb.p[1]-va.p[1];
  double elength = sqrt(dx*dx + dy*dy);
  return fabs(vb.cost-va.cost)/elength;
}    

