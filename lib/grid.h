// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "cell.h"
#include <Eigen/Dense>
#include "map.h"

namespace dsl {

struct EmptyData {};

/**
 * An n-dimenensional grid consisting of abstract "cells", or elements
 * identified by a set of coordinates of type PointT, each cell
 * containing data of type DataT.
 * A grid provides instant access to the elements by maintaining
 * an n-dimensional array of pointers to cells. Cells that are empty,
 * e.g. that are inside obstacles, correspond to null pointers.
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 */

template < class PointT, class DataT = EmptyData >
using Grid = Map< Cell< PointT, DataT >*, PointT::SizeAtCompileTime >;
}
