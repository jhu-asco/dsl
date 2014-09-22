// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this fileSplinecan obtain one at http://mozilla.org/MPL/2.0/.

#include "gridcost3dprm.h"
#include <particleprm.h>
#include <cmath>

using namespace dsl;
#define MAX(a,b) (a>b?a:b)

double GridCost3DPRM::Real(const Vertex &a, const Vertex &b) const
{ 
  const double v = 1.0;
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

  double xs[6] = {x0,y0,z0, v*cos(t0)*sin(p0), v*sin(t0)*sin(p0), v*cos(p0)};
  double xg[6] = {x3,y3,z3, v*cos(t3)*sin(p3), v*sin(t3)*sin(p3), v*cos(p3)};

  System sys(6,3);
  sys.SetBounds(xlb, xub, alb, aub);

  State slb(sys, xlb);
  State sub(sys, xub);
  ParticlePrm prm(slb, sub);
  prm.h = .1;

  Trajectory traj(sys);

  State ss(sys, xs);
  State sg(sys, xg);
  prm.Compute(traj, ss, sg);


  double len = 0;
  for(int i = 0; i < traj.sn; i++)
  {
    
    double dx = traj.states[i+1]->x[0] - traj.states[i]->x[0];
    double dy = traj.states[i+1]->x[1] - traj.states[i]->x[1];
    double dz = traj.states[i+1]->x[2] - traj.states[i]->x[2];
    len += sqrt(dx*dx +dy*dy +dz*dz);
  }
  return len;
}


double GridCost3DPRM::Heur(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = fabs((double)(s1pos[0] - s2pos[0]));
  double dy = fabs((double)(s1pos[1] - s2pos[1]));
  double dz = fabs((double)(s1pos[2] - s2pos[2]));
  return sqrt(dx*dx+dy*dy+dz*dz);//MAX(dz,MAX(dx,dy));
}
