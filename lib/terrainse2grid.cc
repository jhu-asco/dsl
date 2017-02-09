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

using namespace std;
using namespace Eigen;

bool TerrainData::SerializeToString(std::string* str) const{
  int double_size = sizeof(double);
  str->resize(2*double_size);
  str->replace(0,double_size,(char*)&height, double_size);
  str->replace(double_size,double_size,(char*)&traversibility, double_size);
  return true;
}

bool TerrainData::ParseFromString(const std::string& str){
  int n_bytes = sizeof(double);
  std::memcpy(&height,str.c_str(), n_bytes);
  std::memcpy(&traversibility,str.c_str()+n_bytes, n_bytes);
  return true;
}

TerrainData::TerrainData(double height, double traversibility)
:height(height),traversibility(traversibility){

}


TerrainSE2Grid::TerrainSE2Grid(const Map<bool, 3> &cmap, const Map<TerrainData,2>& tmap,
                               const Eigen::Vector3d& cs)
: TerrainSE2GridBase(cmap.xlb(), cmap.xub(), cs,
                     TerrainSE2GridBase::Vectornb(true, false, false),
                     TerrainSE2GridBase::Vectornb(true, true, true)) {

  //Iterate over all cells
  auto fun = [&](int id, const Vector3i& gidx){
    Vector3d cc; //CellCenter
    Vector2i gidx2d = gidx.tail<2>();
    bool gotcenter = CellCenter(gidx, &cc);
    assert(gotcenter);
    bool occ = cmap.Get(cc, false);
    if (!occ) {  //Allocate memory for grid cells if it is not occupied
      this->Set(id,SE2TerrainCell(id, cc, tmap.Get(gidx2d)));
    }
  };
  LoopOver(fun);
}
} //namespace dsl
