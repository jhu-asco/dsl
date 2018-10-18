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

CarGrid::CarGrid(const Map<bool, 3> &cmap,
                 const Vector3d& cs)
    : Grid< Vector3d, Matrix3d >(cmap.xlb, cmap.xub, cs),
      cmap(cmap) {
  for (int k = 0; k < gs[0]; ++k) {
    for (int c = 0; c < gs[1]; ++c) {
      for (int r = 0; r < gs[2]; ++r) {
        // center of cell
        Vector3d x = xlb + Vector3d((k + 0.5) * cs[0], (c + 0.5) * cs[1], (r + 0.5) * cs[2]);
          int id = r*gs[0]*gs[1] + c*gs[0] + k;

        bool occ = cmap.data(x, false);
        if (!occ) {
          values[id] = new SE2Cell(id, x);

          se2_q2g(values[id]->data, values[id]->center);
        }
      }
    }
  }
}


void CarGrid::MakeMap(const Map<bool, 2> &map, Map<bool, 3> &cmap) {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  for (int i = 0; i < cmap.gs[0]; ++i) {
    for (int j = 0; j < cmap.gs[1]; ++j) {
      for (int k = 0; k < cmap.gs[2]; ++k) {
        int id2 = j + cmap.gs[1]*k; //2d index ito map
        assert(id2 < map.nc);
        int id3 = i + cmap.gs[0]*j + cmap.gs[0]*cmap.gs[1]*k; // 3d index into cmap
        assert(id3 < cmap.nc);
        cmap.values[id3] = map.values[id2];
      }
    }
  }
}

void CarGrid::Slice(const Map<bool, 3> &cmap, double a, Map<bool, 2> &map) const {
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
      map.values[id2] = cmap.values[id3];
    }
  }
}


void CarGrid::MakeMap(const CarGeom& geom, const Map<bool, 2> &map, Map<bool, 3> &cmap) {
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

    // dilated map
    //    bool dmap[cmap.gs[1]*cmap.gs[2]];
    //    DilateMap(geom, theta,
    //              cmap.cs[1], cmap.cs[2], cmap.gs[1], cmap.gs[2],
    //              map.values, dmap);
    for (int ix = 0; ix < cmap.gs[1]; ++ix) {
      double x = (ix + 0.5)*cmap.cs[1] + cmap.xlb[1];

      for (int iy = 0; iy < cmap.gs[2]; ++iy) {
        // index into workspace
        int id2 = ix + iy*cmap.gs[1];
        assert(id2 < map.nc);
        // if free continue
        if (!map.values[id2])
          continue;

        double y = (iy + 0.5)*cmap.cs[2] + cmap.xlb[2];

        // index into configuration space
        //        int id3 = ia + ix*cmap.gs[0] + iy*cmap.gs[0]*cmap.gs[1];

        Vector2d p0(x,y); // position of car origin
        for (auto&& dp : points) {
          Vector2d p = p0 + R*dp; // point on the car
          cmap.setData(Vector3d(theta, p[0], p[1]), true);
        }
      }
    }
  }
}

/*
void CarGrid::MakeMap(const CarGeom& geom, const Map<bool, 2> &map, Map<bool, 3>
&cmap) {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  for (int ia = 0; ia < cmap.gs[0]; ++ia) {
    // create a dilated map for a particular angle
    double theta = cmap.xlb[0] + (ia + 0.5) * cmap.cs[0];

    // dilated map
    bool dmap[cmap.gs[1]*cmap.gs[2]];
    DilateMap(geom, theta,
              cmap.cs[1], cmap.cs[2], cmap.gs[1], cmap.gs[2],
              map.values, dmap);
    for (int ix = 0; ix < cmap.gs[1]; ++ix) {
      for (int iy = 0; iy < cmap.gs[2]; ++iy) {
        cmap.values[ia + ix*cmap.gs[0] + iy*cmap.gs[0]*cmap.gs[1]] = dmap[ix +
iy*cmap.gs[1]];
      }
    }
  }
}
*/

 void CarGrid::DilateMap(const CarGeom& geom, double theta,
                         double sx, double sy, int gx, int gy,
                         const bool* data, bool* data_dil) {

   Matrix2x4d verts2d_rotd_pix;
  getRotdVertsInPixWrtOrg(verts2d_rotd_pix, geom.l, geom.b, geom.ox, geom.oy, sx, sy, theta);

  // round of the pixel values of the vertices above such that the rectange
  // formed by the rounded off
  //  vertices surrounds the rotated rectange
  Vector2i org2i_rotd_pix(0,0); // because it's wrt org itself and rounding doesn't matter
  Matrix2x4i verts2i_rotd_pix;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 4; j++)
      verts2i_rotd_pix(i, j) = verts2d_rotd_pix(i, j) > 0 ?
          ceil(verts2d_rotd_pix(i, j)) :
          floor(verts2d_rotd_pix(i, j));

  // Size of kernel is given by the horizontal rectange that bounds the rounded
  // off rectange above
  Vector2i verts2i_rotd_pix_min = verts2i_rotd_pix.rowwise().minCoeff();
  Vector2i verts2i_rotd_pix_max = verts2i_rotd_pix.rowwise().maxCoeff();
  Vector2i size2i_k = verts2i_rotd_pix_max - verts2i_rotd_pix_min;

  // Shift everything such that verts2i_rotd_pix_min is the [0,0] pixel of the
  // kernel
  Matrix2x4i verts2i_rotd_pospix =
      verts2i_rotd_pix.colwise() - verts2i_rotd_pix_min;
  Vector2i org2i_rotd_pospix = org2i_rotd_pix - verts2i_rotd_pix_min;

  // create dilation kernel by filling the inside of the rotated rectanges with
  // zero
  int w_k = size2i_k(0);
  int h_k = size2i_k(1);
  bool data_k[w_k * h_k];
  fillQuad<bool>(data_k,
           size2i_k(0),
           size2i_k(1),
           verts2i_rotd_pospix.cast< double >(),
           1.0);

  // Dilate
  dilate<bool>(data_dil,
         data,
         gx,
         gy,
         data_k,
         size2i_k(0),
         size2i_k(1),
         org2i_rotd_pospix(0),
         org2i_rotd_pospix(1));
 }
}
