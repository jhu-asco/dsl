/*
 * utilsimg.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: subhransu
 */
#include "utilsimg.h"

namespace dsl {

void getRotdVertsInPixWrtOrg(Matrix2x4d& verts2d_rotd_pix, double l,double b, double ox, double oy, double sx, double sy, double theta){
  Vector2d org_m(-ox,-oy);
  Matrix2x4d verts_wrt_center;
  verts_wrt_center.col(0) << -l/2, -b/2;//rb(right back)
  verts_wrt_center.col(1) <<  l/2, -b/2;//rf(right front)
  verts_wrt_center.col(2) <<  l/2,  b/2;//lf(left front)
  verts_wrt_center.col(3) << -l/2,  b/2;//lb(left back)

  //rotate vertices about origin and convert vertices and origin to pixel coordinates
  Transform2d tfm2d= Rotation2Dd(theta)*Translation2d(org_m);
  verts2d_rotd_pix = Scaling(1/sx,1/sy)* (tfm2d*verts_wrt_center);
}

void fillQuad(double* data,int w, int h, Matrix2x4d  verts, double val){

  for(int r=0;r<h;r++){
    for(int c=0;c<w;c++){
      int idx_2d = c + r*w;
      if( inPoly(verts,Vector2d(c,r)))
        data[idx_2d] = val;
      else
        data[idx_2d] = 0;
    }
  }
}

}


