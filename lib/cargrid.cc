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

using namespace dsl;
using namespace Eigen;
using namespace std;

namespace dsl {

using Eigen::Vector3d;
using Eigen::Matrix3d;

CarGrid::CarGrid(const Map<bool, 3> &cmap,
                 const Vector3d& cs) 
    : Grid< Vector3d, Matrix3d >(cmap.xlb, cmap.xub, cs),
      cmap(cmap) {  
  for (int k = 0; k < gs[0]; ++k) {
    for (int c = 0; c < gs[1]; ++c) {
      for (int r = 0; r < gs[2]; ++r) {
        // center of cell
        Vector3d x = xlb + Vector3d((k + 0.5) * cs[0], (c + 0.5) * cs[1], (r + 0.5) * cs[2]);

        bool occ = cmap.Get(x, false);
        if (!occ) {
          int id = r*gs[0]*gs[1] + c*gs[0] + k;
          cells[id] = new SE2Cell(id,
                                  x);
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
        int id3 = i + cmap.gs[0]*j + cmap.gs[0]*cmap.gs[1]*k; // 3d index into cmap        
        cmap.cells[id3] = map.cells[id2];
      }
    }
  }
}


/*
  
  const int& angRes = gs[0];
  for (int k = 0; k < angRes; ++k) {
    // create a dilated map for a particular angle
    double theta = xlb(0) + (k + 0.5) * sa;
    double map_data_dil[map.width * map.height];
    getDilatedMap(map_data_dil, map.data, theta);

    for (int c = 0; c < map.width; ++c) {
      for (int r = 0; r < map.height; ++r) {
        int idx_2d = r * map.width + c; // since data is in row major format
        int idx_3d = r * angRes * map.width + c * angRes +
            k; // 1,2 and 3 dim are a,x and y respectively

        double cost = costScale *
            map_data_dil[idx_2d]; // Cell cost based on angle and geometry of car

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable,
        // so
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[idx_3d] = new SE2Cell(
              xlb + Vector3d((k + 0.5) * sa, (c + 0.5) * sx, (r + 0.5) * sy),
              Vector3d(sa / 2, sx / 2, sy / 2),
              cost);
          se2_q2g(cells[idx_3d]->data, cells[idx_3d]->c);
        }
      }
    }
  }
}

CarGrid::CarGrid(const Map2d &map,
                 double sx,
                 double sy,
                 double sa,
                 double costScale,
                 double maxCost)
  : Grid< 3, Matrix3d >(Vector3d(-M_PI + sa / 2, 0, 0),
                        Vector3d(M_PI + sa / 2, sx * map.width, sy * map.height),
                        Vector3i((int)round(2 * M_PI / sa), map.width, map.height)),
    maxCost(maxCost) {
  const int& angRes = gs[0];
  for (int c = 0; c < map.width; ++c) {
    for (int r = 0; r < map.height; ++r) {
      int idx_2d = r * map.width + c;
      // The pixel value in
      double cost =
          map.data[idx_2d] * costScale; // cell cost = height/occupany/traversability
      for (int k = 0; k < angRes; ++k) {
        int idx_3d = r * angRes * map.width + c * angRes + k;

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable,
        // so
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[idx_3d] = new SE2Cell(
              xlb + Vector3d((k + 0.5) * sa, (c + 0.5) * sx, (r + 0.5) * sy),
              Vector3d(sa / 2, sx / 2, sy / 2),
              cost);
          se2_q2g(cells[idx_3d]->data, cells[idx_3d]->c);
        }
      }
    }
  }
}
*/


void CarGrid::MakeMap(const CarGeom& geom, const Map<bool, 2> &map, Map<bool, 3> &cmap) {

  for (int ia = 0; ia < cmap.gs[0]; ++ia) {
    // create a dilated map for a particular angle
    double theta = cmap.xlb[0] + (ia + 0.5) * cmap.cs[0];

    // dilated map
    bool dmap[cmap.gs[1]*cmap.gs[2]];
    DilateMap(geom, theta,
              cmap.cs[1], cmap.cs[2], cmap.gs[1], cmap.gs[2], 
              map.cells, dmap);
    for (int ix = 0; ix < cmap.gs[1]; ++ix) {
      for (int iy = 0; iy < cmap.gs[2]; ++iy) {
        cmap.cells[ia + ix*cmap.gs[0] + iy*cmap.gs[0]*cmap.gs[1]] = dmap[ix + iy*cmap.gs[1]];
        //cmap.cells[ia + ix*cmap.gs[0] + iy*cmap.gs[0]*cmap.gs[1]] = dmap[iy + ix*cmap.gs[2]];
      }
    }    
  }
}
    

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
