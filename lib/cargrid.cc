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

CarGeom::CarGeom (double l, double b, double ox, double oy) : l(l), b(b), ox(ox), oy(oy) {}

void CarGeom::Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const {
  points.clear();
  for (double x = -l/2; x < l/2; x += cs[0])
    for (double y = -b/2; y < b/2; y += cs[1])
      points.push_back(Eigen::Vector2d(x + ox, y + oy));
}

CarGrid::CarGrid(const Map<bool, 3> &cmap,
                 const Vector3d& cs) 
    : Grid< Vector3d, Matrix3d >(cmap.xlb, cmap.xub, cs,Vector3i(1,0,0)),
      cmap(cmap) {  

  //Allocate memory for grid cells if it is not occupied
  for (int idx_a = 0; idx_a < gs[0]; ++idx_a) {
    for (int idx_x = 0; idx_x < gs[1]; ++idx_x) {
      for (int idx_y = 0; idx_y < gs[2]; ++idx_y) {
        // center of cell
        Vector3i idx(idx_a,idx_x,idx_y);
        Vector3d cc; //CellCenter
        bool gotcenter = CellCenter(cc,idx,true);
        assert(!gotcenter);

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
        cmap.cells[id3] = map.cells[id2];
      }
    }
  }
}

void CarGrid::Slice(const Map<bool, 3> &cmap, double a, Map<bool, 2> &map) {
  map.gs[0] = cmap.gs[1];
  map.gs[1] = cmap.gs[2];
  map.cs[0] = cmap.cs[1];
  map.cs[1] = cmap.cs[2];
  map.cells.resize( round(cmap.nc/cmap.gs[0]));

  int ia = cmap.Index(a, 0);
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


//void CarGrid::MakeMap(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap) {
//  assert(omap.gs[0] == cmap.gs[1]);
//  assert(omap.gs[1] == cmap.gs[2]);
//
//  vector<Vector2d> points;
//  geom.Raster(omap.cs/2, points);
//
//  Matrix2d R;
//
//  int dim_a(0); //dimension for angle
//  int dim_x(1); //dimension for angle
//  int dim_y(2); //dimension for angle
//
//  for (int idx_a = 0; idx_a < cmap.gs[dim_a]; ++idx_a) {
//    // dilate map for a particular angle
//    double theta = cmap.CellCenterIth(idx_a,dim_a);
//
//    // make a rotation matrix
//    double ct = cos(theta);
//    double st = sin(theta);
//    R(0,0) = ct; R(0,1) = -st;
//    R(1,0) = st; R(1,1) = ct;
//
//    for (int idx_x = 0; idx_x < cmap.gs[dim_x]; ++idx_x) {
//      double x = cmap.CellCenterIth(idx_x,dim_x);
//      for (int idx_y = 0; idx_y < cmap.gs[2]; ++idx_y) {
//        // index into workspace
//        int id_omap = omap.Id(Vector2i(idx_x,idx_y));
//        assert(id_omap < omap.nc);
//
//        if (!omap.cells[id_omap])
//          continue;        // if free continue
//
//        double y = cmap.CellCenterIth(idx_y,dim_y);
//
//        Vector2d p0(x,y); // position of car origin
//        for (auto&& dp : points) {
//          Vector2d p = p0 + R*dp; // point on the car
//          cmap.Set(Vector3d(theta, p[0], p[1]), true);
//        }
//      }
//    }
//  }
//}


void CarGrid::MakeMap(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap) {
  assert(map.gs[0] == cmap.gs[1]);
  assert(map.gs[1] == cmap.gs[2]);

  int dim_a = 0; //The dimension corresponding to angle

  for (int idx_a = 0; idx_a < cmap.gs[dim_a]; ++idx_a) {
    // create a dilated map for a particular angle
    double theta = cmap.CellCenterIth(idx_a,dim_a);

    // dilated map
    vector<bool> dmap(omap.nc);
    DilateMap(geom, theta, cmap.cs[1], cmap.cs[2], cmap.gs[1], cmap.gs[2],omap.cells, dmap);

    for (int idx_x = 0; idx_x < cmap.gs[1]; ++idx_x) {
      for (int idx_y = 0; idx_y < cmap.gs[2]; ++idx_y) {
        int id_omap = omap.Id( Vector2i(idx_x, idx_y) );
        int id_cmap = cmap.Id( Vector3i(idx_a, idx_x, idx_y) );
        cmap.cells[id_cmap] = dmap[id_omap];
      }
    }
  }

}

 void CarGrid::DilateMap(const CarGeom& geom, double theta,
                         double sx, double sy, int gx, int gy, 
                         const vector<bool>& data, vector<bool>& data_dil) {
   
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
  vector<bool> data_k(w_k * h_k);
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
