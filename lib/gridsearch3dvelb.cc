// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <math.h>
#include <stdlib.h>
#include "gridsearch3dvelb.h"

using namespace std;
using namespace dsl;

#define _USE_MATH_DEFINES
#define MAX(a,b) (a>b?a:b)

void GridSearch3DVelB::GetTrajectory(const Vertex &from, const Vertex &to, GridPath3DPlusTime &path) const
{
  const double rayTraceStep = 0.05;
  int numSteps = 1.0/rayTraceStep;
  double len = 0;
  path.pos = (double*) realloc(path.pos, numSteps*3*sizeof(double));
  path.times = (double*) realloc(path.times, numSteps*sizeof(double));
  path.count = numSteps;

  const double smoothness = 0.8;
  int* s1pos = (int*)from.data;
  int* s2pos = (int*)to.data;
  double x0 = s1pos[0];
  double y0 = s1pos[1];
  double z0 = s1pos[2];
  double t0 = s1pos[3]*M_PI/180.;
  double p0 = s1pos[4]*M_PI/180.;
  double x3 = s2pos[0];
  double y3 = s2pos[1];
  double z3 = s2pos[2];
  double t3 = s2pos[3]*M_PI/180.;
  double p3 = s2pos[4]*M_PI/180.;
  double x1,y1,z1,x2,y2,z2;


  x1 = x0 + smoothness*cos(t0)*sin(p0);
  y1 = y0 + smoothness*sin(t0)*sin(p0);
  z1 = z0 + smoothness*cos(p0);

  x2 = x3 - smoothness*cos(t3)*sin(p3);
  y2 = y3 - smoothness*sin(t3)*sin(p3);
  z2 = z3 - smoothness*cos(p3);

  //cout << "v0: " << x0 << " " << y0 << " " << z0
  //<< " v1: " << x1 << " " << y1 << " " << z1
  //<< " v2: " << x2 << " " << y2 << " " << z2
  //<< " v3: " << x3 << " " << y3 << " " << z3 << endl;

  double xi,yi,zi,xil,yil,zil;
  xil = x0;
  yil = y0;
  zil = z0;
  for(int j = 0; j < numSteps; j++)
  {
    double t = j*rayTraceStep;
    xi = (1-t)*(1-t)*(1-t)*x0 + 3*(1-t)*(1-t)*t*x1 + 3*(1-t)*t*t*x2 + t*t*t*x3;
    yi = (1-t)*(1-t)*(1-t)*y0 + 3*(1-t)*(1-t)*t*y1 + 3*(1-t)*t*t*y2 + t*t*t*y3;
    zi = (1-t)*(1-t)*(1-t)*z0 + 3*(1-t)*(1-t)*t*z1 + 3*(1-t)*t*t*z2 + t*t*t*z3;
    if(int(xi) < 0 || int(xi) >= length || int(yi) < 0 || int(yi) >= width || int(zi) < 0 || int(zi) >= height || (int)GetCost((int)xi,(int)yi,(int)zi) == DSL3D_OCCUPIED)
    {
      path.len = -1;
      return;
    }
    path.times[j] = t;
    path.pos[j*3] = xi;
    path.pos[j*3+1] = yi;
    path.pos[j*3+2] = zi;
    double dx = xi-xil;
    double dy = yi-yil;
    double dz = zi-zil;
    len += sqrt(dx*dx+dy*dy+dz*dz);
    xil = xi;
    yil = yi;
    zil = zi;
  }
  path.len = len;
  return;
}
