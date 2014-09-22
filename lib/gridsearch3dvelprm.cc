// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <map>
#include "gridsearch3dvelprm.h"
#include "particleprm.h"

using namespace std;
using namespace dsl;
#define _USE_MATH_DEFINES
#define MAX(a,b) ((a)>(b)?(a):(b))

map< PathInfo, Trajectory* > pathCache;

void GridSearch3DVelPRM::GetTrajectory(const Vertex &from, const Vertex &to, GridPath3DPlusTime &path) const
{
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

  PathInfo pi;
  pi.dx = x3-x0;
  pi.dy = y3-y0;
  pi.dz = z3-z0;
  pi.t1 = t0;
  pi.t2 = t3;
  pi.p1 = p0;
  pi.p2 = p3;

  std::map< PathInfo, Trajectory* >::iterator it;
  
  Trajectory *cachedPath;
  if((it = pathCache.find(pi)) != pathCache.end())
  {
    //cout << "Cached path found" << endl;
    cachedPath = it->second;
  }
  else
  {
    //cout << "Cached path NOT found" << endl;
    double xs[6] = {0,0,0, this->v*cos(t0)*sin(p0), this->v*sin(t0)*sin(p0), this->v*cos(p0)};
    double xg[6] = {pi.dx,pi.dy,pi.dz, this->v*cos(t3)*sin(p3), this->v*sin(t3)*sin(p3), this->v*cos(p3)};

    System sys(6,3);
    sys.SetBounds(xlb, xub, alb, aub);
 
    State slb(sys, xlb);
    State sub(sys, xub);
    ParticlePrm prm(slb, sub);
    prm.h = .5;
 
    Trajectory *traj = new Trajectory(sys);
 
    State ss(sys, xs);
    State sg(sys, xg);
    if(!prm.Compute(*traj, ss, sg))
    {
      traj->sn = 0;
      traj->states = NULL;
    }
    pathCache[pi] = traj;
    cachedPath = traj; 
  } 

  if(cachedPath->states == NULL)
  {
    path.len = -1;
    return;
  }

  //cout << "traj.sn: " << traj.sn << endl; 
  path.pos = (double*) realloc(path.pos, (cachedPath->sn+1)*3*sizeof(double));
  path.times = (double*) realloc(path.times, (cachedPath->sn+1)*sizeof(double));
  path.count = cachedPath->sn;

  double len = 0;
  for(int i = 0; i < cachedPath->sn; i++)
  {
    double x = cachedPath->states[i]->x[0] + x0;
    double y = cachedPath->states[i]->x[1] + y0;
    double z = cachedPath->states[i]->x[2] + z0;
    if((int)x < 0 || (int)y < 0 || (int) z < 0 || (int) x >= length || (int)y >= width || (int)z >= height || GetCost((int)x, (int)y, (int)z) == DSL3D_OCCUPIED)
    {
      path.len = -1;
      return;
    }
    path.times[i] = cachedPath->states[i]->t;
    path.pos[i*3] = x;
    path.pos[i*3+1] = y;
    path.pos[i*3+2] = z;
    double dx = cachedPath->states[i+1]->x[0] - cachedPath->states[i]->x[0];
    double dy = cachedPath->states[i+1]->x[1] - cachedPath->states[i]->x[1];
    double dz = cachedPath->states[i+1]->x[2] - cachedPath->states[i]->x[2];
    len += sqrt(dx*dx +dy*dy +dz*dz);
  }

  path.len = len;
}
