// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "cost.h"

using namespace dsl;

double Cost::Real(const Vertex &va, const Vertex &vb) const
{
  return Heur(va, vb) + 1e-10;
}
