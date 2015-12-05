// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gridcost3d.h"
#include <cmath>

using namespace dsl;

#define MAX(a,b) (a>b?a:b)

double GridCost3D::Real(const Cell3d &a, const Cell3d &b) const
{
  double dx = a.p[0] - b.p[0];
  double dy = a.p[1] - b.p[1];
  double dz = a.p[2] - b.p[2];
  return sqrt(dx*dx + dy*dy + dz*dz)*(1+ 0.5*a.cost + 0.5*b.cost);
}


double GridCost3D::Heur(const Cell3d &a, const Cell3d &b) const
{
  double dx = fabs((double)(a.p[0] - b.p[0]));
  double dy = fabs((double)(a.p[1] - b.p[1]));
  double dz = fabs((double)(a.p[2] - b.p[2]));
  return MAX(dz,MAX(dx,dy));
}
