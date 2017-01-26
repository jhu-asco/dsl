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
using namespace Eigen;
using namespace std;

//SE2CellGrid::Vectornb(true,false,false) says which dimensions are wrapped
CarGrid::CarGrid(const Map<bool, 3> &cmap, const Vector3d& cs)
    : SE2CellGrid(cmap.xlb(), cmap.xub(), cs, SE2CellGrid::Vectornb(true,false,false)){

  //Iterate over all cells
  auto fun = [&](int id, const Vector3i& gidx){
    Vector3d cc; //CellCenter
    bool gotcenter = CellCenter(gidx, &cc);
    assert(gotcenter);
    bool occ = cmap.Get(cc, false);
    if (!occ) {  //Allocate memory for grid cells if it is not occupied
      Matrix3d g;
      se2_q2g(g, cc);
      //this->set_cells(id,SE2Cell(id, cc, g));
      //cells_[id].reset(new SE2Cell(id, cc));
    }
  };
  LoopOver(fun);
}
} /* namespace dsl */
