// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "cargrid.h"
#include "utils.h"

#include <iostream>
#include "utilsimg.h"
#include <thread>

namespace dsl {

using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;
using std::vector;

//SE2CellGrid::Vectornb(true,false,false) says which dimensions are wrapped
CarGrid::CarGrid(const Map<bool, 3> &cmap, const Vector3d& cs)
    : SE2CellGrid(cmap.xlb, cmap.xub, cs, SE2CellGrid::Vectornb(true,false,false)){

  //Allocate memory for grid cells if it is not occupied
  for (int idx_a = 0; idx_a < gs[0]; ++idx_a) {
    for (int idx_x = 0; idx_x < gs[1]; ++idx_x) {
      for (int idx_y = 0; idx_y < gs[2]; ++idx_y) {
        // center of cell
        Vector3i idx(idx_a,idx_x,idx_y);
        Vector3d cc; //CellCenter
        bool gotcenter = CellCenter(cc,idx,true);
        assert(gotcenter);

        bool occ = cmap.Get(cc, false);
        if (!occ) {
          int id = Id(idx);
          cells[id].reset(new SE2Cell(id, cc));
          se2_q2g(cells[id]->data, cells[id]->c);
        }
      }
    }
  }
}
}
