// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "search.h"
#include "fibheap.h"

namespace dsl
{
  // these are needed by fibheap
  int FIBHEAPKEY_SIZE = 2*sizeof(double);
  fibheapkey_t FIBHEAPKEY_MIN = (void*)(double[2]){-INF, -INF};

  extern "C" int fibkey_compare(fibheapkey_t a, fibheapkey_t b)
  {
    assert(a); assert(b);
    double* af = (double *)a;
    double* bf = (double *)b;
    if ((af[0] < bf[0]) || (af[0] == bf[0] && af[1] < bf[1]))
      return -1;
    if (af[0] == bf[0] && af[1] == bf[1])
      return 0;
    return 1;
  }
}
