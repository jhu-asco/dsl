// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this fileSplinecan obtain one at http://mozilla.org/MPL/2.0/.

#include "gridcost3dspline.h"
#include "spline.h"
#include <cmath>

using namespace dsl;

#define MAX(a,b) (a>b?a:b)

double GridCost3DSpline::Real(const Vertex &a, const Vertex &b) const
{ 
  const double timeStep = 0.05;
  const double smoothness = 0.1;
  const double v = 1.0;
  double totalTime = 0;
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
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

  std::vector<double> times;
  std::vector<double> pointsx;
  std::vector<double> pointsy;
  std::vector<double> pointsz;
  
  x1 = x0 + smoothness*cos(t0)*sin(p0);
  y1 = y0 + smoothness*sin(t0)*sin(p0);
  z1 = z0 + smoothness*cos(p0);

  x2 = x3 - smoothness*cos(t3)*sin(p3);
  y2 = y3 - smoothness*sin(t3)*sin(p3);
  z2 = z3 - smoothness*cos(p3);

  pointsx.push_back(x0);
  pointsy.push_back(y0);
  pointsz.push_back(z0);
  pointsx.push_back(x1);
  pointsy.push_back(y1);
  pointsz.push_back(z1);
  pointsx.push_back(x2);
  pointsy.push_back(y2);
  pointsz.push_back(z2);
  pointsx.push_back(x3);
  pointsy.push_back(y3);
  pointsz.push_back(z3);
  times.push_back(0);
  for(int i = 1; i < pointsx.size(); i++)
  {
    double dx = pointsx[i] - pointsx[i-1];
    double dy = pointsy[i] - pointsy[i-1];
    double dz = pointsz[i] - pointsz[i-1];
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    totalTime += dist/v;
    times.push_back(totalTime);
  }
  /* Create the spline interpolating the position over time */
  Spline<double, double> spx(times, pointsx);
  Spline<double, double> spy(times, pointsy);
  Spline<double, double> spz(times, pointsz);

  double len = 0;
  for(int i = 1; i < totalTime/timeStep; i++)
  {
    double dx = spx[i*timeStep] - spx[(i-1)*timeStep];
    double dy = spy[i*timeStep] - spy[(i-1)*timeStep];
    double dz = spz[i*timeStep] - spz[(i-1)*timeStep];
    len += sqrt(dx*dx +dy*dy +dz*dz);
  }

  return len;
}


double GridCost3DSpline::Heur(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = fabs((double)(s1pos[0] - s2pos[0]));
  double dy = fabs((double)(s1pos[1] - s2pos[1]));
  double dz = fabs((double)(s1pos[2] - s2pos[2]));
  return sqrt(dx*dx+dy*dy+dz*dz);//MAX(dz,MAX(dx,dy));
}
