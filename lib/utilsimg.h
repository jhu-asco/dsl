// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra subhransu.kumar.mishra@gmail.com
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_UTILSIMG_H_
#define DSL_LIB_UTILSIMG_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
namespace dsl {

using Transform2d = Eigen::Transform< double, 2, Eigen::Affine >;
using Matrix2x4d = Eigen::Matrix< double, 2, 4 >;
using Matrix2x4i = Eigen::Matrix< int, 2, 4 >;

/**
 * Takes in the dimensions of a rectangle, position of its origin and a rotation
 * angle theta and returns the
 * position of the 4 corners(right-back, right-front, left-front and left-back)
 * of the rotated rectangle
 * with respect to the origin in pixel units.
 * @param verts2d_rotd_pix coordinates of the 4 corners(rb,rf,lf,lb) in order
 * @param l Dimension along the x axis of the car
 * @param b Dimension along the y axis of the car
 * @param ox x coordinate origin of the car with respect to the rectangle center
 * @param oy y coordinate origin of the car with respect to the rectangle center
 * @param sx Pixel width in meters in x direction(horizontal)
 * @param sy Pixel width in meters in y direction(vertical)
 * @param theta Angle by which the rectangle is rotated
 */
void getRotdVertsInPixWrtOrg(Matrix2x4d& verts2d_rotd_pix,
                             double l,
                             double b,
                             double ox,
                             double oy,
                             double sx,
                             double sy,
                             double theta);

/**
 * Scales up an image with the same scaling factor in x and y direction
 * @param map_scaled pointer to scaled image. It needs to be allocated before
 * @param map pointer to the original image
 * @param w Width of the map
 * @param h Height of the map
 * @param scale Scaling factor
 */
template < typename T >
void scaleMap(T* map_scaled, const T* map, int w, int h, int scale) {
  int ws = scale * w;
  int hs = scale * h;

  for (int rs = 0; rs < hs; rs++) {
    for (int cs = 0; cs < ws; cs++) {
      int c = round(cs / scale);
      int r = round(rs / scale);

      int idx = c + r * w;
      int idxs = cs + rs * ws;

      map_scaled[idxs] = map[idx];
    }
  }
}

/**
 * Scales up an image with the same scaling factor in x and y direction
 * @param map_scaled Reference to scaled map
 * @param map Original map
 * @param w Width of the map
 * @param h Height of the map
 * @param scale Scaling factor
 */
template < typename T >
void scaleMap(std::vector<T>& map_scaled, const std::vector<T>& map, int w, int h, int scale) {
  int ws = scale * w;
  int hs = scale * h;

  map_scaled.resize(ws*hs);
  for (int rs = 0; rs < hs; rs++) {
    for (int cs = 0; cs < ws; cs++) {
      int c = round(cs / scale);
      int r = round(rs / scale);

      int idx = c + r * w;
      int idxs = cs + rs * ws;

      map_scaled[idxs] = map[idx];
    }
  }
}

/**
 * Checks if a given point is in a polygon
 * @param verts A Matrix of vertices of the polygon. Each column represents a
 * vertex.
 * @param pt The point to checked
 * @return
 */
template < typename T, int m, int n >
int inPoly(Eigen::Matrix< T, m, n > verts, Eigen::Matrix< T, m, 1 > pt) {
  int i, j, c = 0;
  for (i = 0, j = n - 1; i < n; j = i++) {
    if (((verts(1, i) > pt(1)) != (verts(1, j) > pt(1))) &&
        (pt(0) < (verts(0, j) - verts(0, i)) * (pt(1) - verts(1, i)) /
                 (verts(1, j) - verts(1, i)) +
             verts(0, i)))
      c = !c;
  }
  return c;
}

/**
 * Takes in an image and end points of a line and adds that line of certain
 * pixel width
 * @param map The input image
 * @param w The width of the image
 * @param h The height of the image
 * @param p1 End point1 for the line
 * @param p2 End point2 for the line
 * @param lval the pixel value for the line
 * @param lw linewidth in pixels
 */
template < typename T >
void addLine( T* map, int w, int h, Eigen::Vector2d p1, Eigen::Vector2d p2, T lval, double lw) {
  Eigen::Vector2d n = Eigen::Rotation2Dd(M_PI / 2) * (p2 - p1).normalized();
  Matrix2x4d verts2d;
  verts2d.col(0) = p1 + n * lw / 2;
  verts2d.col(1) = p1 - n * lw / 2;
  verts2d.col(2) = p2 - n * lw / 2;
  verts2d.col(3) = p2 + n * lw / 2;

  for (int r = 0; r < h; r++) {
    for (int c = 0; c < w; c++) {
      int idx = c + r * w;
      if (inPoly(verts2d, Eigen::Vector2d(c, r)))
        map[idx] = lval;
    }
  }
}

/**
 * Takes in an image and end points of a line and adds that line of certain
 * pixel width
 * @param map The input image as a vector
 * @param w The width of the image
 * @param h The height of the image
 * @param p1 End point1 for the line in pixel coordinates
 * @param p2 End point2 for the line in pixel coordinates
 * @param lval the pixel value for the line
 * @param lw linewidth in pixels
 */
template < typename T >
void addLine( std::vector<T>& map, int w, int h, Eigen::Vector2d p1, Eigen::Vector2d p2, T lval, double lw) {
  Eigen::Vector2d n = Eigen::Rotation2Dd(M_PI / 2) * (p2 - p1).normalized();
  Matrix2x4d verts2d;
  verts2d.col(0) = p1 + n * lw / 2;
  verts2d.col(1) = p1 - n * lw / 2;
  verts2d.col(2) = p2 - n * lw / 2;
  verts2d.col(3) = p2 + n * lw / 2;

  for (int r = 0; r < h; r++) {
    for (int c = 0; c < w; c++) {
      int idx = c + r * w;
      if (inPoly(verts2d, Eigen::Vector2d(c, r)))
        map[idx] = lval;
    }
  }
}

/**
 * Takes in an input image, a kernel image a and the origin of kernel image and
 * produces a dilated image.
 * @param data_dil Dilated image
 * @param data Input image
 * @param w Width of the input image
 * @param h Height of the input image
 * @param data_k Kernel image
 * @param w_k Width of the kernel
 * @param h_k Height of the kernel
 * @param ox_k x coordinate of the origin of the kernel image
 * @param oy_k y coordinate of the origin of the kernel image
 */
template < typename T >
void dilate(std::vector<T>& data_dil,
            const std::vector<T>&  data,
            int w,
            int h,
            const std::vector<T>& data_k,
            int w_k,
            int h_k,
            int ox_k,
            int oy_k) {
  std::vector< T > prod(h_k * w_k);

  // Visit each pixel in the main image
  for (int r = 0; r < h; r++) {
    int r_rel = r - oy_k;
    for (int c = 0; c < w; c++) {
      int c_rel = c - ox_k;
      int id = c + r * w;
      // for each pixel in the main image lay the kernel on the original image
      //  such that the origin of kernel(only 1s and 0s) coincides with pixel in
      //  question
      //  Then for pixel in question visit all the elements in the kernel. Each
      //  element of kernel image
      //  has a corresponding pixel in main image(except at the borders).
      //  Multiply those pair of values
      //  together. Take the max of all those values and that becomes the pixel
      //  value of the dilated image.
      //  What this does is checks if the any element of kernel is overlaid over
      //  an obstacle.
      for (int r_k = 0; r_k < h_k; r_k++) {
        for (int c_k = 0; c_k < w_k; c_k++) {
          int id_k = c_k + r_k * w_k;
          int id_roi = c_k + c_rel + (r_k + r_rel) * w;
          if (c_k + c_rel >= 0 && c_k + c_rel < w && r_k + r_rel >= 0 &&
              r_k + r_rel < h)
            prod[id_k] = data_k[id_k] * data[id_roi];
          else
            prod[id_k] = 0;
        }
      }
      data_dil[id] = *(std::max_element< typename std::vector< T >::iterator >(
          prod.begin(), prod.end()));

    }
  }
}


/**
 * Fills a quadrilateral with a given value
 * @param data pointer to the img data
 * @param w Width of the img
 * @param h Height of the img
 * @param verts The 4 vertices of the quadrilateral
 * @param val the value that is to be fill in. Rest is 0
 */
template < typename T = bool>
    void fillQuad(std::vector<T>& data, int w, int h, Matrix2x4d verts, T val) {
  for (int r = 0; r < h; r++) {
    for (int c = 0; c < w; c++) {
      int idx_2d = c + r * w;
      if (inPoly(verts, Eigen::Vector2d(c, r)))
        data[idx_2d] = val;
      else
        data[idx_2d] = 0;
    }
  }
}

}
#endif /* DSL_UTILSIMG_H */
