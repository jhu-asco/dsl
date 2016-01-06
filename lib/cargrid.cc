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

CarGrid::CarGrid(double l,double b, double ox, double oy,
                 int width, int height, const double* map, double sx, double sy, double sa,
                 double costScale,double maxCost)
                 :Grid<3, Matrix3d>(Vector3d(-M_PI + sa/2, 0, 0),
                                    Vector3d(M_PI + sa/2, sx*width, sy*height),
                                    Vector3i((int)round(2*M_PI/sa), width, height))
                 ,maxCost(maxCost),l_(l),b_(b),ox_(ox),oy_(oy){

  const int &angRes = gs[0];
  for (int k = 0; k < angRes; ++k){
    //create a dilated map for a particular angle
    double theta = xlb(0) + (k + 0.5)*sa;
    double map_dil[width*height];
    getDilatedMap(map_dil,map,theta);

    for (int c = 0; c < width; ++c){
      for (int r = 0; r < height; ++r){
        int idx_2d = r*width + c;// since data is in row major format
        int idx_3d = r*angRes*width + c*angRes + k; //1,2 and 3 dim are a,x and y respectively

        double cost = costScale * map_dil[idx_2d];//Cell cost based on angle and geometry of car

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so
        // they shouldn't be added to the list of cells
        if (cost < maxCost){
          cells[idx_3d] = new SE2Cell(xlb + Vector3d((k + 0.5)*sa, (c + 0.5)*sx, (r + 0.5)*sy),
                                      Vector3d(sa/2, sx/2, sy/2), cost);
          se2_q2g(cells[idx_3d]->data, cells[idx_3d]->c);
        }
      }
    }
  }

}

CarGrid::CarGrid(int width, int height, const double *map,
                 double sx, double sy, double sa, double costScale,
                 double maxCost) :
          Grid<3, Matrix3d>(Vector3d(-M_PI+ sa/2, 0, 0), Vector3d(M_PI+ sa/2, sx*width, sy*height),
                            Vector3i((int)round(2*M_PI/sa), width, height)),
                            maxCost(maxCost),l_(0),b_(0),ox_(0),oy_(0) {
  const int &angRes = gs[0];
  for (int c = 0; c < width; ++c){
    for (int r = 0; r < height; ++r){
      int idx_2d = r*width + c;
      //The pixel value in
      double cost = map[idx_2d]*costScale; // cell cost = height/occupany/traversability
      for (int k = 0; k < angRes; ++k){
        int idx_3d = r*angRes*width + c*angRes + k;

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so
        // they shouldn't be added to the list of cells
        if (cost < maxCost)
        {
          cells[idx_3d] = new SE2Cell(xlb + Vector3d((k + 0.5)*sa, (c + 0.5)*sx, (r + 0.5)*sy),
                                  Vector3d(sa/2, sx/2, sy/2), cost);
          se2_q2g(cells[idx_3d]->data, cells[idx_3d]->c);
        }
      }
    }
  }

}



void
CarGrid::getDilatedMap(double* data_dil, const double* data, double theta){
  double sx=cs(1), sy=cs(2);
  Matrix2x4d verts2d_rotd_pix;
  getRotdVertsInPixWrtOrg(verts2d_rotd_pix, l_, b_, ox_, oy_, sx, sy,theta);

  //round of the pixel values of the vertices above such that the rectange formed by the rounded off
  //  vertices surrounds the rotated rectange
  Vector2i org2i_rotd_pix; org2i_rotd_pix.setZero();//because it's wrt org itself and rounding doesn't matter
  Matrix2x4i verts2i_rotd_pix;
  for(int i=0;i<2;i++)
    for(int j=0;j<4;j++)
      verts2i_rotd_pix(i,j) = verts2d_rotd_pix(i,j)>0?ceil(verts2d_rotd_pix(i,j)):floor(verts2d_rotd_pix(i,j));

  //Size of kernel is given by the horizontal rectange that bounds the rounded off rectange above
  Vector2i verts2i_rotd_pix_min = verts2i_rotd_pix.rowwise().minCoeff();
  Vector2i verts2i_rotd_pix_max = verts2i_rotd_pix.rowwise().maxCoeff();
  Vector2i size2i_k = verts2i_rotd_pix_max - verts2i_rotd_pix_min;

  //Shift everything such that verts2i_rotd_pix_min is the [0,0] pixel of the kernel
  Matrix2x4i verts2i_rotd_pospix = verts2i_rotd_pix.colwise() - verts2i_rotd_pix_min;
  Vector2i   org2i_rotd_pospix   = org2i_rotd_pix - verts2i_rotd_pix_min;

  //create dilation kernel by filling the inside of the rotated rectanges with zero
  int w_k=size2i_k(0);
  int h_k=size2i_k(1);
  double data_k[w_k*h_k];
  fillQuad(data_k,size2i_k(0), size2i_k(1),verts2i_rotd_pospix.cast<double>(),1.0);

  //Dilate
  dilate(data_dil, data, gs[1], gs[2], data_k, size2i_k(0), size2i_k(1),
         org2i_rotd_pospix(0), org2i_rotd_pospix(1));
}


}
