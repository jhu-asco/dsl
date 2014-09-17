// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "gridsearch3d.h"
#include "spline.h"

using namespace std;
using namespace dsl;

#define MAX(a,b) (a>b?a:b)

// cell neighbor connectivity
// this can be either 6 (faster,less accurate) or 26 (slower, more accurate) 
#define NBR_COUNT 26
#define SQRT2 1.414213562373095
#define SQRT3 1.73205080757

// the order of these matters
static int NBR_OFFSETS[NBR_COUNT*3] = {0,0,-1, -1,0,-1, -1,-1,-1, 0,-1,-1, 1,-1,-1, 1,0,-1, 1,1,-1, 0,1,-1, -1,1,-1, -1,0,0, -1,-1,0, 0,-1,0, 1,-1,0, 1,0,0, 1,1,0, 0,1,0, -1,1,0, 0,0,1, -1,0,1, -1,-1,1, 0,-1,1, 1,-1,1, 1,0,1, 1,1,1, 0,1,1, -1,1,1};
static double NBR_COSTS_[NBR_COUNT] = {1.0, SQRT2, SQRT3, SQRT2, SQRT3, SQRT2, SQRT3, SQRT2, SQRT3, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, SQRT3, SQRT2, SQRT3, SQRT2, SQRT3, SQRT2, SQRT3};


GridPath3D::GridPath3D() : pos(0), count(0), len(0) { }

GridPath3D::GridPath3D(const GridPath3D& cp)
{
  if(cp.pos)
  {
    pos = (double*)malloc(3*cp.count*sizeof(double));
    memcpy(pos, cp.pos, 3*cp.count*sizeof(double));
  }
  else
  {
    pos = 0;
  }
  len = cp.len;
  count = cp.count;
}

GridPath3D::~GridPath3D() 
{
  if(pos)
  {
    free(pos);
    pos = 0;
  }
}

GridPath3DPlusTime::GridPath3DPlusTime() : GridPath3D(), times(0) {}
GridPath3DPlusTime::GridPath3DPlusTime(const GridPath3DPlusTime& cp)
 : GridPath3D(cp)
{
  if(cp.times)
  {
    times = (double*)malloc(cp.count*sizeof(double));
    memcpy(times, cp.times, cp.count*sizeof(double));
  }
  else
  {
    times = 0;
  }
}

GridPath3DPlusTime::~GridPath3DPlusTime()
{
  if(times)
    free(times);
}


void GridPath3DPlusTime::AppendPath(GridPath3DPlusTime &p)
{
  if(p.times)
  {
    double* cur_times = this->times;
    this->times = (double*) malloc((this->count + p.count)*sizeof(double));

    for(int i = 0; i < this->count; i++)
    {
      this->times[i] = cur_times[i];
    }
    for(int i = 0; i < p.count; i++)
    {
      this->times[this->count+i] = p.times[i] + this->times[this->count-1];
    }
    if(this->len > 0 && cur_times)
      free(cur_times);
  }
  this->GridPath3D::AppendPath(p);
}

void GridPath3D::AppendPath(GridPath3D &p)
{
  double* cur_path = this->pos;
  this->pos = (double*) malloc((this->count + p.count)*3*sizeof(double));

  for(int i = 0; i < 3*this->count; i++)
  {
    this->pos[i] = cur_path[i];
  }
  for(int i = 0; i < 3*p.count; i++)
  {
    this->pos[3*this->count+i] = p.pos[i];
  }
  if(this->len > 0 && cur_path)
    free(cur_path);

  this->len += p.len;
  this->count += p.count;  
}

GridSearch3D::GridSearch3D(int length, int width, int height, const double *map, double scale) : 
  Search(graph, cost),
  length(length),
  width(width),
  height(height), 
  scale(scale)
{
  int x, y, z, i, nbr, nx, ny, nz, ni;
  int* pos;

  int s = length*width*height;

  // initialize map
  this->map = new double[s];
  if (map)
    memcpy(this->map, map, s*sizeof(double));
  else 
    memset(this->map, 0, s*sizeof(double));
  
  //cout << "Grid map initialized" << endl;
  // this will be used to map grid coordinates to dsl graph vertices
  vertexMap = new Vertex*[s];
  
  // create all vertices
  i = 0;
  for(z = 0; z < height; ++z){
    for (y = 0; y < width; ++y) {
      for (x = 0; x < length; ++x, ++i) {
        pos = new int[3];
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
        vertexMap[i] = new Vertex(pos);
        graph.AddVertex(*vertexMap[i]);
      }            //      this->vertices.push_back();
    }
  }

  //cout << "Created vertices" << endl;

  // create all edges
  i = 0;
  for(z = 0; z < height; ++z){
    for (y = 0; y < width; ++y) {
      for (x = 0; x < length; ++x, ++i) {    
        for (nbr = 0; nbr < NBR_COUNT; ++nbr) {
	  nx = x + NBR_OFFSETS[3*nbr];
	  ny = y + NBR_OFFSETS[3*nbr+1];
	  nz = z + NBR_OFFSETS[3*nbr+2];
	  if (nx < 0 || ny < 0 || nz < 0 || nx >= length || ny >= width || nz >= height)
	    continue;
	  ni = nz*width*length + ny*length + nx;
          // create an edge b/n vertices at pos. i and ni

          Vertex *from = vertexMap[i];
          Vertex *to = vertexMap[ni];
        
        

          Edge* edge = new Edge(from, to, scale*(MAX(this->map[i], this->map[ni]) + NBR_COSTS_[nbr]));
        

          graph.AddEdge(*edge);

          //        from->out[edge->id] = edge;
          //        to->in[edge->id] = edge;
          //        edges.push_back(edge);
        }
      }
    }
    //cout << z << endl;
  }
}


void GridSearch3D::AddEdge(int x1, int y1, int z1, int x2, int y2, int z2)
{  
  Vertex *from = GetVertex(x1, y1, z1);
  Vertex *to = GetVertex(x2, y2, z2);
  if (!from || !to)
    return;
  
  double dx = x2-x1;
  double dy = y2-y1;
  double dz = z2-z1;
  Edge* edge = new Edge(from, to, sqrt(dx*dx + dy*dy + dz*dz)); 
  graph.AddEdge(*edge);
}


GridSearch3D::~GridSearch3D()
{
  for (int i = 0; i < length*width*height; ++i) {
    delete[] (int*)vertexMap[i]->data;
    vertexMap[i]->data = 0;
  }
  delete[] vertexMap;
  delete[] map;
}


double GridSearch3D::GetCost(int x, int y, int z) const
{
  if(x < 0 || y < 0 || z < 0 || x >= length || y >= width || z >= height)
    return DSL3D_OCCUPIED;

  return map[z*length*width + y*length + x];
}


Vertex* GridSearch3D::GetVertex(int x, int y, int z) const
{
  if (x < 0 || x >= length || y < 0 || y >= width || z < 0 || z >= height)
    return 0;

  return vertexMap[z*length*width + y*length + x];
}

void GridSearch3D::RemoveVertex(int x, int y, int z)
{
  Vertex *v = GetVertex(x, y, z);
  if (v)
    graph.RemoveVertex(*v);
}


void GridSearch3D::SetCost(int x, int y, int z, double cost)
{
  std::map<int, Edge*>::iterator ein, eout;

  assert(cost >= 0);
  
  int i = z*length*width + y*length + x;
 
  // if the cost has not changed simply return
  if (Eq(cost, this->map[i]))
    return;
  this->map[i] = cost;
  
  // fix all connected edges
  ein = vertexMap[i]->in.begin();
  eout = vertexMap[i]->out.begin();
  for (;ein != vertexMap[i]->in.end(); ein++)
    ChangeCost(*ein->second, cost > 10 ? 10000 : 0);

  for (;eout != vertexMap[i]->out.end(); eout++)
    ChangeCost(*eout->second, cost > 10 ? 10000 : 0);
}

void GridSearch3D::SetMap(const double *map)
{
  int z, x, y, i;
  i = 0;
  for(z = 0; z < height; ++z){
    for (y = 0; y < width; ++y) {
      for (x = 0; x < length; ++x, ++i) {
        SetCost(x, y, z, map[i]);
      }
    }
  }
  memcpy(this->map, map, length*width*height*sizeof(double));
}


void GridSearch3D::SetStart(int x, int y, int z)
{
  Search::SetStart(*(GetVertex(x,y,z)));
}

void GridSearch3D::SetGoal(int x, int y, int z)
{
  Search::SetGoal(*(GetVertex(x,y,z)));
}

void GridSearch3D::Plan(GridPath3D& path)
{
  int count;
  double len = 0;
  Vertex *cur = start;
  int *pos1 = 0, *pos0 = 0;
  int i;

  count = Search::Plan();
  path.pos = (double*)realloc(path.pos, count*3*sizeof(double));
  path.count = count;
  for (i = 0; i < count; ++i) {
    pos1 = (int*)cur->data;
    path.pos[3*i] = (double)pos1[0];   
    path.pos[3*i+1] = (double)pos1[1];   
    path.pos[3*i+2] = (double)pos1[2];   
    //memcpy(&path.pos[3*i], pos1, 3*sizeof(double));
    if (i > 0) {
      len += sqrt((pos1[0]-pos0[0])*(pos1[0]-pos0[0]) + (pos1[1]-pos0[1])*(pos1[1]-pos0[1]) + (pos1[2]-pos0[2])*(pos1[2]-pos0[2]));
    }
    pos0 = pos1;
    cur = cur->next;
  }
  path.len = len;
}

#define RAY_TRACE_STEP 0.05

void GridSearch3D::SmoothPathOptCost(const GridPath3D &path, GridPath3DPlusTime &smoothPath, double v, double timeStep) const
{
  srand(time(NULL));
  double totalTimeOrig = 0;
  std::vector<double> times;
  std::vector<double> pointsx;
  std::vector<double> pointsy;
  std::vector<double> pointsz;
  for(int i = 0; i < path.count; i++)
  {
    pointsx.push_back(path.pos[3*i]);
    pointsy.push_back(path.pos[3*i+1]);
    pointsz.push_back(path.pos[3*i+2]);
    if(i>0)
    {
      double dx = path.pos[3*i] - path.pos[3*(i-1)];    
      double dy = path.pos[3*i+1] - path.pos[3*(i-1)+1];
      double dz = path.pos[3*i+2] - path.pos[3*(i-1)+2];
      double dist = sqrt(dx*dx + dy*dy + dz*dz);
      //cout << path.pos[3*i] << " "
      // << path.pos[3*i+1] << " "
      // << path.pos[3*i+2] << " "
      // << totalTime << endl;
      totalTimeOrig += dist/v;
    }
    times.push_back(totalTimeOrig);
  }

  /* Create the spline interpolating the position over time */
  Spline<double, double> spxOrig(times, pointsx);
  Spline<double, double> spyOrig(times, pointsy);
  Spline<double, double> spzOrig(times, pointsz);

  Spline<double, double> spxOpt;
  Spline<double, double> spyOpt;
  Spline<double, double> spzOpt;
  double totalTimeOpt;
  double minCost = 999999999;
  for(int i = 0; i < 200000; i++)
  {
    std::vector<double> timesTest;
    std::vector<double> pointsxTest;
    std::vector<double> pointsyTest;
    std::vector<double> pointszTest;

    // Always add first point
    timesTest.push_back(0);
    pointsxTest.push_back(path.pos[0]);
    pointsyTest.push_back(path.pos[1]);
    pointszTest.push_back(path.pos[2]);

    // Randomly choose points along path to add as spline ctrl points
    int curPt = 1;
    int numCtrlPts = 0;
    while(curPt < path.count-2)
    {
      int pt = (rand() % (path.count - curPt)) + curPt;
        
      double dx = path.pos[3*pt] - pointsxTest.at(pointsxTest.size()-1);    
      double dy = path.pos[3*pt+1] - pointsyTest.at(pointsyTest.size()-1);    
      double dz = path.pos[3*pt+2] - pointszTest.at(pointszTest.size()-1);    
      double dist = sqrt(dx*dx + dy*dy + dz*dz);
 
      timesTest.push_back(timesTest.at(timesTest.size()-1) + dist/v);
      pointsxTest.push_back(path.pos[3*pt]);
      pointsyTest.push_back(path.pos[3*pt+1]);
      pointszTest.push_back(path.pos[3*pt+2]);
      
      numCtrlPts++;
      curPt = pt+1;
    }
    if(numCtrlPts < 1)
      continue;

    // Always add last point
    double dx = path.pos[3*(path.count-1)] - pointsxTest.at(pointsxTest.size()-1);    
    double dy = path.pos[3*(path.count-1)+1] - pointsyTest.at(pointsyTest.size()-1);    
    double dz = path.pos[3*(path.count-1)+2] - pointszTest.at(pointszTest.size()-1);    
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
 
    timesTest.push_back(timesTest.at(timesTest.size()-1) + dist/v);
    pointsxTest.push_back(path.pos[3*(path.count-1)]);
    pointsyTest.push_back(path.pos[3*(path.count-1)+1]);
    pointszTest.push_back(path.pos[3*(path.count-1)+2]);

    // Calculate path cost
    Spline<double, double> spxTest(timesTest, pointsxTest);
    Spline<double, double> spyTest(timesTest, pointsyTest);
    Spline<double, double> spzTest(timesTest, pointszTest);
    
    double ctrlPtCost = 5*(numCtrlPts+2.)/path.count;
    double pathDistCost = 0;
    double tScale = timesTest.at(timesTest.size()-1)/totalTimeOrig;

    for(double i = 0; i < totalTimeOrig; i += timeStep)
    {
      double x = spxTest[tScale*i];
      double y = spyTest[tScale*i];
      double z = spzTest[tScale*i];
      double dx = spxOrig[i] - x;
      double dy = spyOrig[i] - y;
      double dz = spzOrig[i] - z;

      // Make sure path doesn't collide
      if ( GetCost((int)x, (int)y, (int)z) == DSL3D_OCCUPIED)
      { 
        pathDistCost = -10.; 
        break;
      }

      pathDistCost += sqrt(dx*dx + dy*dy + dz*dz);
    }
    pathDistCost /= (totalTimeOrig/timeStep);

    if(pathDistCost > 0 && ctrlPtCost + pathDistCost < minCost)
    {
      cout << "NumCtrlPts: " << numCtrlPts <<  " CtrlPtCost: " << ctrlPtCost  << " PathDistCost: " << pathDistCost << " Total Cost: " << (ctrlPtCost + pathDistCost) << endl;     
      minCost = ctrlPtCost + pathDistCost;
      spxOpt = spxTest;
      spyOpt = spyTest;
      spzOpt = spzTest;
      totalTimeOpt = timesTest.at(timesTest.size()-1);
    }    
  }
  
  // Convert to dsl path
  int count = totalTimeOpt/timeStep + 1;
  smoothPath.pos = (double*) realloc(smoothPath.pos, count*3*sizeof(double));
  smoothPath.times = (double*) realloc(smoothPath.times, count*sizeof(double));
  smoothPath.count = count;
  for(int i = 0; i < count-1; i++)
  {
    smoothPath.pos[3*i] = spxOpt[i*timeStep];
    smoothPath.pos[3*i+1] = spyOpt[i*timeStep];
    smoothPath.pos[3*i+2] = spzOpt[i*timeStep];
    smoothPath.times[i] = i*timeStep;
  }

  smoothPath.pos[3*(count-1)] = spxOpt[totalTimeOpt];
  smoothPath.pos[3*(count-1)+1] = spyOpt[totalTimeOpt];
  smoothPath.pos[3*(count-1)+2] = spzOpt[totalTimeOpt];
  smoothPath.times[count-1] = totalTimeOpt;

  double len = 0;
  for(int i = 0; i < count-1; i++)
  {
    double dx = smoothPath.pos[(i+1)*3] - smoothPath.pos[i*3];
    double dy = smoothPath.pos[(i+1)*3+1] - smoothPath.pos[i*3+1];
    double dz = smoothPath.pos[(i+1)*3+2] - smoothPath.pos[i*3+2];
    len += sqrt(dx*dx + dy*dy + dz*dz);
  }
  smoothPath.len = len;

}

void GridSearch3D::SmoothPathSpline(const GridPath3D &path, GridPath3DPlusTime &smoothPath, double v, double timeStep) const
{
  double totalTime = 0;
  std::vector<double> times;
  std::vector<double> pointsx;
  std::vector<double> pointsy;
  std::vector<double> pointsz;
  for(int i = 0; i < path.count; i++)
  {
    pointsx.push_back(path.pos[3*i]);
    pointsy.push_back(path.pos[3*i+1]);
    pointsz.push_back(path.pos[3*i+2]);
    if(i>0)
    {
      double dx = path.pos[3*i] - path.pos[3*(i-1)];      
      double dy = path.pos[3*i+1] - path.pos[3*(i-1)+1];      
      double dz = path.pos[3*i+2] - path.pos[3*(i-1)+2];
      double dist = sqrt(dx*dx + dy*dy + dz*dz);    
      totalTime += dist/v;  
    }
    times.push_back(totalTime);
  }

  /* Create the spline interpolating the position over time */
  Spline<double, double> spx(times, pointsx);
  Spline<double, double> spy(times, pointsy);
  Spline<double, double> spz(times, pointsz);

  int count = times.at(times.size()-1)/timeStep + 1;
  smoothPath.pos = (double*) realloc(smoothPath.pos, count*3*sizeof(double));
  smoothPath.times = (double*) realloc(smoothPath.times, count*sizeof(double));
  smoothPath.count = count;
  for(int i = 0; i < count-1; i++)
  {
    //cout << spx[i*timeStep] << " " 
    //<< spy[i*timeStep] << " " 
    //<< spz[i*timeStep] << " " 
    //<< i*timeStep << endl;

    smoothPath.pos[3*i] = spx[i*timeStep];
    smoothPath.pos[3*i+1] = spy[i*timeStep];
    smoothPath.pos[3*i+2] = spz[i*timeStep];
    smoothPath.times[i] = i*timeStep;
  }
 
  smoothPath.pos[3*(count-1)] = spx[times.at(times.size()-1)];
  smoothPath.pos[3*(count-1)+1] = spy[times.at(times.size()-1)];
  smoothPath.pos[3*(count-1)+2] = spz[times.at(times.size()-1)];
  smoothPath.times[count-1] = times.at(times.size()-1);
  
  double len = 0;
  for(int i = 0; i < count-1; i++)
  {
    double dx = smoothPath.pos[(i+1)*3] - smoothPath.pos[i*3];
    double dy = smoothPath.pos[(i+1)*3+1] - smoothPath.pos[i*3+1];
    double dz = smoothPath.pos[(i+1)*3+2] - smoothPath.pos[i*3+2];
    len += sqrt(dx*dx + dy*dy + dz*dz);
  }
  smoothPath.len = len;
}



void GridSearch3D::SmoothPathBezier(const GridPath3D &path, GridPath3D &smoothPath, double smoothness) const
{
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, n;
  int numSteps = 1.0/RAY_TRACE_STEP;
  double count =  (numSteps*(path.count-2)+2);
  double len = 0;
  smoothPath.pos = (double*) realloc(smoothPath.pos, count*3*sizeof(double));
  smoothPath.count = count;

  // Interpolate with cubic bezier curve
  for(int i = 0; i < path.count-2; i++)
  {
    x0 = path.pos[3*i];
    y0 = path.pos[3*i+1];
    z0 = path.pos[3*i+2];
    x3 = path.pos[3*(i+1)];
    y3 = path.pos[3*(i+1)+1];
    z3 = path.pos[3*(i+1)+2];
    
    x1 = x3 - x0;
    y1 = y3 - y0;
    z1 = z3 - z0;
    n = sqrt(x1*x1 + y1*y1 + z1*z1);
    x1 *= smoothness/n;    
    y1 *= smoothness/n;    
    z1 *= smoothness/n;
    x1 += x0;   
    y1 += y0;   
    z1 += z0;   
    
    x2 = path.pos[3*(i+2)] - x3;
    y2 = path.pos[3*(i+2)+1] - y3;
    z2 = path.pos[3*(i+2)+2] - z3;
    n = sqrt(x2*x2 + y2*y2 + z2*z2);
    x2 *= -smoothness/n;    
    y2 *= -smoothness/n;    
    z2 *= -smoothness/n;
    x2 += x3;   
    y2 += y3;   
    z2 += z3;   
    for(int j = 0; j < numSteps; j++)
    {
      double t = j*RAY_TRACE_STEP;
      int idx = (i)*numSteps + j;
      smoothPath.pos[idx*3] = (1-t)*(1-t)*(1-t)*x0 + 3*(1-t)*(1-t)*t*x1 + 3*(1-t)*t*t*x2 + t*t*t*x3;
      smoothPath.pos[idx*3+1] = (1-t)*(1-t)*(1-t)*y0 + 3*(1-t)*(1-t)*t*y1 + 3*(1-t)*t*t*y2 + t*t*t*y3;
      smoothPath.pos[idx*3+2] = (1-t)*(1-t)*(1-t)*z0 + 3*(1-t)*(1-t)*t*z1 + 3*(1-t)*t*t*z2 + t*t*t*z3;
    }
  }
  
  for(int i = -6; i < 0; i ++)
  {
    smoothPath.pos[smoothPath.count*3+i] = path.pos[path.count*3+i];
  } 

  for(int i = 0; i < count-1; i++)
  {
    double dx = smoothPath.pos[(i+1)*3] - smoothPath.pos[i*3];
    double dy = smoothPath.pos[(i+1)*3+1] - smoothPath.pos[i*3+1];
    double dz = smoothPath.pos[(i+1)*3+2] - smoothPath.pos[i*3+2];
    len += sqrt(dx*dx + dy*dy + dz*dz);
  }
  smoothPath.len = len;
  //std::cout << "Count: " << count << " Path count: " << path.count << std::endl;

}


void GridSearch3D::OptPath(const GridPath3D &path, GridPath3D &optPath) const
{
  double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, n, d;
  double dx0, dy0, dz0;
  double dx1, dy1, dz1;
  int i;
  int count = 1;
  double len = 0;
  double pos[3*path.count];

  if (path.len == 2) {
    optPath.pos = (double*)realloc(optPath.pos, 6*sizeof(double));
    memcpy(optPath.pos, path.pos, 6*sizeof(double));
    optPath.count = 2;
    optPath.len = path.len;
    return;
  }
  
  x2 = 0;
  y2 = 0;
  z2 = 0;

  x0 = path.pos[0];
  y0 = path.pos[1];
  z0 = path.pos[2];
  x1 = path.pos[3];
  y1 = path.pos[4];
  z1 = path.pos[5];
  dx0 = x1 - x0;
  dy0 = y1 - y0;
  dz0 = z1 - z0;
  n = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
  dx0 /= n;
  dy0 /= n;
  dz0 /= n;

  pos[0] = x0;
  pos[1] = y0;
  pos[2] = z0;
  
  for (i = 1; i < path.count - 1; ++i) {
    x1 = path.pos[3*i];
    y1 = path.pos[3*i + 1];
    z1 = path.pos[3*i + 2];
    x2 = path.pos[3*(i+1)];
    y2 = path.pos[3*(i+1) + 1];
    z2 = path.pos[3*(i+1) + 2];
    dx1 = x2 - x0;
    dy1 = y2 - y0;
    dz1 = z2 - z0;
    n = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
    assert(n > .7);
    dx1 /= n;
    dy1 /= n;
    dz1 /= n;
    if (fabs((dx0-dx1)*(dx0-dx1) + (dy0-dy1)*(dy0-dy1) + (dz0-dz1)*(dz0-dz1)) > .00001) {
      n = sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2) + (z0-z2)*(z0-z2));
      for (d = RAY_TRACE_STEP; d < n; d += RAY_TRACE_STEP) {
	x = dx1*d;
	y = dy1*d;
        z = dz1*d;
	if (((int)map[((int)(z0 + z))*length*width + ((int)(y0 + y))*length + (int)(x0 + x)]) == DSL3D_OCCUPIED) {
	  pos[3*count] = x1;
	  pos[3*count + 1] = y1;
	  pos[3*count + 2] = z1;
	  count++;
	  len += sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
	  x0 = x1;
	  y0 = y1;
          z0 = z1;
	  break;
	}
      }
      dx0 = x2 - x0;
      dy0 = y2 - y0;
      dz0 = z2 - z0;
      n = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
      dx0 /= n;
      dy0 /= n;
      dz0 /= n;
    }
  }
  
  pos[3*count] = x2;
  pos[3*count + 1] = y2;
  pos[3*count + 2] = z2;
  len += sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));
  count++;
  optPath.pos = (double*)realloc(optPath.pos, count*3*sizeof(double));
  memcpy(optPath.pos, pos, count*3*sizeof(double));
  optPath.count = count;
  optPath.len = len;
}

