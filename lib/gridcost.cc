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


double GridCost::Real(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = s1pos[0] - s2pos[0];
  double dy = s1pos[1] - s2pos[1];
  return sqrt(dx*dx + dy*dy);
}


double GridCost::Heur(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = fabs((double)(s1pos[0] - s2pos[0]));
  double dy = fabs((double)(s1pos[1] - s2pos[1]));
  if (dx > dy)
    return dx;
  else
    return dy;
}
