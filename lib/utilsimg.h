/*
 * utilsimg.h
 *
 *  Created on: Jan 4, 2016
 *      Author: subhransu
 */

#ifndef DSL_UTILSIMG_H
#define DSL_UTILSIMG_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
namespace dsl {

using namespace Eigen;
using namespace std;

typedef Transform<double,2,Affine> Transform2d;
typedef Matrix<double,2,4> Matrix2x4d;
typedef Matrix<int,2,4> Matrix2x4i;


void getRotdVertsInPixWrtOrg(Matrix2x4d& verts2d_rotd_pix,
                             double l,double b, double ox, double oy, double sx, double sy, double theta);

void fillQuad(double* data,int w, int h, Matrix2x4d  verts, double val);


template <typename T>
void scaleMap(T* map_scaled, T* map, int w, int h, int scale){
  int ws = scale*w;
  int hs = scale*h;

  for( int rs=0;rs<hs;rs++){
    for(int cs=0;cs<ws;cs++){
      int c = round(cs/scale);
      int r = round(rs/scale);

      int idx = c+r*w;
      int idxs = cs+rs*ws;

      map_scaled[idxs] = map[idx];
    }
  }
}

template< typename T, int m, int n>
int inPoly(Matrix<T,m,n> verts, Matrix<T,m,1> pt){
  int i, j, c = 0;
  for (i = 0, j = n-1; i < n; j = i++) {
    if ( ((verts(1,i)>pt(1)) != (verts(1,j)>pt(1))) &&
        (pt(0) < (verts(0,j)-verts(0,i)) * (pt(1)-verts(1,i)) / (verts(1,j)-verts(1,i)) + verts(0,i)) )
      c = !c;
  }
  return c;
}


template <typename T>
void addLine(T* map, int w, int h, Vector2d p1, Vector2d p2, T lval, double lw){
  Vector2d n = Rotation2Dd(M_PI/2) * (p2-p1).normalized();
  Matrix2x4d verts2d;
  verts2d.col(0) = p1 + n*lw/2;
  verts2d.col(1) = p1 - n*lw/2;
  verts2d.col(2) = p2 - n*lw/2;
  verts2d.col(3) = p2 + n*lw/2;

  for( int r=0;r<h;r++){
    for(int c=0;c<w;c++){
      int idx = c+r*w;
      if(inPoly(verts2d, Vector2d(c,r)))
        map[idx] = lval;
    }
  }
}

template<typename T>
void dilate(T* data_dil, T* data, int w, int h, T* data_k, int w_k, int h_k, int ox_k, int oy_k){
  vector<double> prod(h_k*w_k);

  //Visit each pixel in the main image
  for(int r=0; r<h; r++){
    int r_rel = r - oy_k;
    for(int c=0; c <w; c++){
      int c_rel = c - ox_k;
      int id=c+r*w;
      //for each pixel in the main image lay the kernel on the original image
      //  such that the origin of kernel(only 1s and 0s) coincides with pixel in question
      //  Then for pixel in question visit all the elements in the kernel. Each element of kernel image
      //  has a corresponding pixel in main image(except at the borders). Multiply those pair of values
      //  together. Take the max of all those values and that becomes the pixel value of the dilated image.
      //  What this does is checks if the any element of kernel is overlaid over an obstacle.
      for(int r_k=0; r_k < h_k; r_k++){
        for(int c_k=0; c_k < w_k; c_k++){
          int id_k=c_k+r_k*w_k;
          int id_roi = c_k+c_rel + (r_k+r_rel)*w;
          if(c_k + c_rel>=0 && c_k + c_rel<w && r_k+r_rel>=0 && r_k+r_rel<h  )
            prod[id_k] = data_k[id_k] * data[id_roi];
          else
            prod[id_k] = 0;
        }
      }
      data_dil[id] = *(max_element<typename vector<T>::iterator> (prod.begin(),prod.end()));
    }
  }
}

}
#endif /* DSL_UTILSIMG_H */
