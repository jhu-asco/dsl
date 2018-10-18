// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "cell.h"
#include "indexed_array.h"

namespace dsl {

struct EmptyData {};

template < class PointT, class DataT = EmptyData >
using Grid = IndexedArray< PointT, Cell< PointT, DataT >* >;
}
