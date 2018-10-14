// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra subhransu.kumar.mishra@gmail.com
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "utilsimg.h"

namespace dsl {

void getRotdVertsInPixWrtOrg(Matrix2x4d& verts2d_rotd_pix,
                             double l,
                             double b,
                             double ox,
                             double oy,
                             double sx,
                             double sy,
                             double theta) {
  Vector2d org_m(-ox, -oy);
  Matrix2x4d verts_wrt_centr;
  verts_wrt_centr.col(0) << -l / 2, -b / 2; // rb(right back)
  verts_wrt_centr.col(1) << l / 2, -b / 2;  // rf(right front)
  verts_wrt_centr.col(2) << l / 2, b / 2;   // lf(left front)
  verts_wrt_centr.col(3) << -l / 2, b / 2;  // lb(left back)

  // rotate vertices about origin and convert vertices and origin to pixel
  // coordinates
  Transform2d tfm2d = Rotation2Dd(theta) * Translation2d(org_m);
  verts2d_rotd_pix = Scaling(1 / sx, 1 / sy) * (tfm2d * verts_wrt_centr);
}

}
