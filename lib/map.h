// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "grid.h"
#include <iostream>
#include <fstream>

namespace dsl {

template < class T, int n >
using Map = IndexedArray< Eigen::Matrix< double, n, 1 >, T >;

template < typename T, int n >
static void save(const Map< T, n >& map, const char* filename) {
  std::ofstream fs(filename, std::fstream::out | std::ios::binary);
  assert(fs.is_open());

  for (int i = 0; i < map.xlb.size(); ++i)
    fs.write((char*)&map.xlb[i], sizeof(double));
  for (int i = 0; i < map.xub.size(); ++i)
    fs.write((char*)&map.xub[i], sizeof(double));
  for (int i = 0; i < map.gs.size(); ++i)
    fs.write((char*)&map.gs[i], sizeof(double));

  fs.write((char*)map.values, map.nc * sizeof(T));
  fs.close();
}

template < typename T, int n >
static Map< T, n >* load(const char* filename) {
  using Vectornd = Eigen::Matrix< double, n, 1 >;
  using Vectorni = Eigen::Matrix< int, n, 1 > ;

  std::ifstream fs(filename, std::fstream::in | std::ios::binary);
  assert(fs.is_open());
  Vectornd xlb, xub;
  Vectorni gs;

  for (int i = 0; i < xlb.size(); ++i)
    fs.read((char*)&xlb[i], sizeof(double));
  for (int i = 0; i < xub.size(); ++i)
    fs.read((char*)&xub[i], sizeof(double));
  for (int i = 0; i < gs.size(); ++i)
    fs.read((char*)&gs[i], sizeof(double));

  Map< T, n >* map = new Map< T, n >(xlb, xub, gs);
  fs.read((char*)map->values, map->nc * sizeof(T));
  fs.close();
  return map;
}
}
