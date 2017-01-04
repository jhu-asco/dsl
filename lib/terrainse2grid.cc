// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <terrainse2grid.h>
#include "utils.h"

#include <iostream>
#include "utilsimg.h"
#include <thread>

namespace dsl {

using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;
using std::vector;

TerrainData::TerrainData(double height, double traversibility)
:height(height),traversibility(traversibility){

}


TerrainSE2Grid::TerrainSE2Grid(const Map<bool, 3> &cmap, const Map<TerrainData,2>& tmap,
                               const Eigen::Vector3d& cs)
: TerrainSE2GridBase(cmap.xlb, cmap.xub, cs,
                     TerrainSE2GridBase::Vectornb(true, false, false),
                     TerrainSE2GridBase::Vectornb(true, true, true)) {

  //Iterate over all cells
  auto fun = [&](int id, const Vector3i& gidx){
    Vector3d cc; //CellCenter
    Vector2i gidx2d = gidx.tail<2>();
    bool gotcenter = CellCenter(cc,gidx);assert(gotcenter);
    bool occ = cmap.Get(cc, false);
    if (!occ) {  //Allocate memory for grid cells if it is not occupied
      cells[id].reset(new TerrainCell(id, cc));
      cells[id]->data = tmap.Get(gidx2d);
    }
  };
  LoopOver(fun);

}
} //namespace dsl
