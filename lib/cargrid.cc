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

namespace dsl {

using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;

using std::vector;

CarGrid::CarGrid(const Map3b &cmap,
                 const Vector3d& cs)
    : SE2Grid(cmap.xlb, cmap.xub, cs),
      cmap(cmap) {
  for (int k = 0; k < gs[0]; ++k) {
    for (int c = 0; c < gs[1]; ++c) {
      for (int r = 0; r < gs[2]; ++r) {
        // center of cell
        Vector3d x = xlb + Vector3d((k + 0.5) * cs[0], (c + 0.5) * cs[1], (r + 0.5) * cs[2]);
          int id = r*gs[0]*gs[1] + c*gs[0] + k;

        bool occ = cmap.data(x, false);
        if (!occ) {
          cells[id] = new SE2Cell(id, x);

          se2_q2g(cells[id]->data, cells[id]->center);
        }
      }
    }
  }
}


void CarGrid::MakeMap(const Map2b &map, Map3b &cmap) {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  for (int i = 0; i < cmap.gs[0]; ++i) {
    for (int j = 0; j < cmap.gs[1]; ++j) {
      for (int k = 0; k < cmap.gs[2]; ++k) {
        int id2 = j + cmap.gs[1]*k; //2d index ito map
        assert(id2 < map.nc);
        int id3 = i + cmap.gs[0]*j + cmap.gs[0]*cmap.gs[1]*k; // 3d index into cmap
        assert(id3 < cmap.nc);
        cmap.cells[id3] = map.cells[id2];
      }
    }
  }
}

void CarGrid::Slice(const Map3b &cmap, double a, Map2b &map) const {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  int ia = cmap.index(a, 0);
  for (int ix = 0; ix < cmap.gs[1]; ++ix) {
    for (int iy = 0; iy < cmap.gs[2]; ++iy) {
      // index into workspace
      int id2 = ix + iy*cmap.gs[1];
      assert(id2 < map.nc);
        // index into configuration space
      int id3 = ia + ix*cmap.gs[0] + iy*cmap.gs[0]*cmap.gs[1];
      assert(id3 < cmap.nc);
      map.cells[id2] = cmap.cells[id3];
    }
  }
}


void CarGrid::MakeMap(const CarGeom& geom, const Map2b &map, Map3b &cmap) {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  vector<Vector2d> points;
  geom.raster(map.cs, points);

  Matrix2d R;

  for (int ia = 0; ia < cmap.gs[0]; ++ia) {
    // dilate map for a particular angle
    double theta = (ia + 0.5) * cmap.cs[0] + cmap.xlb[0];

    // make a rotation matrix
    double ct = cos(theta);
    double st = sin(theta);
    R(0,0) = ct; R(0,1) = -st;
    R(1,0) = st; R(1,1) = ct;

    for (int ix = 0; ix < cmap.gs[1]; ++ix) {
      double x = (ix + 0.5)*cmap.cs[1] + cmap.xlb[1];

      for (int iy = 0; iy < cmap.gs[2]; ++iy) {
        // index into workspace
        int id2 = ix + iy*cmap.gs[1];
        assert(id2 < map.nc);
        // if free continue
        if (!map.cells[id2])
          continue;

        double y = (iy + 0.5)*cmap.cs[2] + cmap.xlb[2];

        Vector2d p0(x,y); // position of car origin
        for (auto&& dp : points) {
          Vector2d p = p0 + R*dp; // point on the car
          cmap.setData(Vector3d(theta, p[0], p[1]), true);
        }
      }
    }
  }
}

}
