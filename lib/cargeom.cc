// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "cargeom.h"

namespace dsl {

CarGeom::CarGeom (double l, double b, double ox, double oy, double sb){
  set(l,b,ox,oy,sb);
}

CarGeom::CarGeom (const Vector5d& lboxoysb){
  set(lboxoysb);
}

void CarGeom::set(double l, double b, double ox, double oy, double sb){
  l_ = l;
  b_ = b;
  ox_ = ox;
  oy_ =oy;
  sb_ = sb;

  le_ = l_ + 2*sb;
  be_ = b_ + 2*sb;
}
void CarGeom::set(const Vector5d& lboxoysb){
  l_ = lboxoysb[0];
  b_ = lboxoysb[1];
  ox_ = lboxoysb[2];
  oy_ = lboxoysb[3];
  sb_ = lboxoysb[4];

  le_ = l_ + 2*sb_;
  be_ = b_ + 2*sb_;
}

void CarGeom::Raster(const Eigen::Vector2d &cs, std::vector<Eigen::Vector2d> &points) const {
  points.clear();
  for (double x = -le_/2; x < le_/2; x += cs[0])
    for (double y = -be_/2; y < be_/2; y += cs[1])
      points.push_back(Eigen::Vector2d(x + ox_, y + oy_));
}

double CarGeom::l() const {return l_;}
double CarGeom::b() const {return b_;}
double CarGeom::ox() const {return ox_;}
double CarGeom::oy() const {return oy_;}
double CarGeom::sb() const {return sb_;}
double CarGeom::le() const {return le_;}
double CarGeom::be() const {return be_;}

void CarGeom::GetTrueCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const{

  Matrix2x4d verts_in_org;
  verts_in_org.col(0) << -l_/2 -ox_, -b_/2 -oy_; // rb(right back)
  verts_in_org.col(1) <<  l_/2 -ox_, -b_/2 -oy_;  // rf(right front)
  verts_in_org.col(2) <<  l_/2 -ox_,  b_/2 -oy_;   // lf(left front)
  verts_in_org.col(3) << -l_/2 -ox_,  b_/2 -oy_;  // lb(left back)

  Transform2d rotation; rotation = Rotation2Dd(theta);
  Matrix2x4d verts_rotd =  rotation * verts_in_org;

  vertices.resize(4);
  for(int i=0; i<4; i++)
    vertices[i] = verts_rotd.col(i);
}

void CarGeom::GetSafeCorners(std::vector<Eigen::Vector2d>& vertices, double theta) const{

  Matrix2x4d verts_in_org;
  verts_in_org.col(0) << -le_/2 -ox_, -be_/2 -oy_; // rb(right back)
  verts_in_org.col(1) <<  le_/2 -ox_, -be_/2 -oy_;  // rf(right front)
  verts_in_org.col(2) <<  le_/2 -ox_,  be_/2 -oy_;   // lf(left front)
  verts_in_org.col(3) << -le_/2 -ox_,  be_/2 -oy_;  // lb(left back)

  Transform2d rotation; rotation = Rotation2Dd(theta);
  Matrix2x4d verts_rotd =  rotation * verts_in_org;

  vertices.resize(4);
  for(int i=0; i<4; i++)
    vertices[i] = verts_rotd.col(i);
}

} //namespace dsl
