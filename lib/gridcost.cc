// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gridcost.h"
#include <cmath>

using namespace dsl;


double GridCost::Real(const Cell2d &a, const Cell2d &b) const
{
  int dx = a.p[0] - b.p[0];
  int dy = a.p[1] - b.p[1];
  return sqrt(dx*dx + dy*dy)*(1 + 0.5*a.cost + 0.5*b.cost);
}


double GridCost::Heur(const Cell2d &a, const Cell2d &b) const
{
  double dx = fabs((double)(a.p[0] - b.p[0]));
  double dy = fabs((double)(a.p[1] - b.p[1]));
  if (dx > dy)
    return dx;
  else
    return dy;
}
