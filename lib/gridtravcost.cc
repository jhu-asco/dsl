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
#include "gridsearch.h"

using namespace std;
using namespace dsl;


double GridTravCost::Real(const Vertex &va, const Vertex &vb) const
{
  VertexGridData* adat = ( VertexGridData* )va.data;
  VertexGridData* bdat = ( VertexGridData* )vb.data;
  int dx = bdat->p[0]-adat->p[0];
  int dy = bdat->p[1]-adat->p[1];
  double elength = sqrt(dx*dx + dy*dy);
  return fabs(bdat->cost-adat->cost)/elength;
}    

