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
#include "travedgecost.h"

using namespace std;
using namespace dsl;


double TravEdgeCost::CalcEdgeCost(double v1cost, double v2cost, double elength)
{
  assert(elength > 0);
  return fabs(v2cost- v1cost)/elength;
}

