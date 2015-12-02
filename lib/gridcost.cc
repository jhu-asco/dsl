// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gridcost.h"
#include "gridsearch.h"
#include <cmath>

using namespace dsl;


double GridCost::Real(const Vertex &a, const Vertex &b) const
{
  VertexGridData* va_dat = (VertexGridData*)a.data;
  VertexGridData* vb_dat = (VertexGridData*)b.data;
  double vcosta = va_dat->cost;
  double vcostb = vb_dat->cost;
  int dx = va_dat->p[0] - vb_dat->p[0];
  int dy = va_dat->p[1] - vb_dat->p[1];
  return sqrt(dx*dx + dy*dy)*(1+0.5*vcosta+0.5*vcostb);
}


double GridCost::Heur(const Vertex &a, const Vertex &b) const
{
  VertexGridData* va_dat = (VertexGridData*)a.data;
  VertexGridData* vb_dat = (VertexGridData*)b.data;
  int* s1pos = va_dat->p;
  int* s2pos = vb_dat->p;
  double dx = fabs((double)(s1pos[0] - s2pos[0]));
  double dy = fabs((double)(s1pos[1] - s2pos[1]));
  if (dx > dy)
    return dx;
  else
    return dy;
}
