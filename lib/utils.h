// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sys/time.h>
#include <Eigen/Dense>
#include <vector>
#include "grid.h"

namespace dsl {

Map2b fromPPM(const char* filename, const Eigen::Vector2d& cs);

void toPPM(const Map2b& map,
               const char* filename,
               const std::vector< Eigen::Vector2d >* path = 0);

void toPPM(const char* map, int width, int height, const char* filename);

char* fromPPM(int &width, int &height, const char* filename);

void timer_start(struct timeval* time);

long timer_us(struct timeval* time);

double fangle(double a);

void se2_q2g(Eigen::Matrix3d& m, const Eigen::Vector3d& q);

void se2_g2q(Eigen::Vector3d& q, const Eigen::Matrix3d& m);

void se2_inv(Eigen::Matrix3d& mi, const Eigen::Matrix3d& m);

void se2_exp(Eigen::Matrix3d& m, const Eigen::Vector3d& v, double tol = 1e-16);

void replaceExt(std::string& s, const std::string& newExt);

/**
 * Function that returns a sign for object of any class as long as the operator
 * - and operator < are defined
 * for that class
 * @param val the object whose sign we need to check
 * @return sign of the object
 */
template < typename T >
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
}
