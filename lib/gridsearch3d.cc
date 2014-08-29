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

//static int NBR_OFFSETS[8*2] = {1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1};
//static double NBR_COSTS_[8] = {1.0, 1.0, 1.0, 1.0, SQRT2, SQRT2, SQRT2, SQRT2};


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
  this->map = const_cast<double*>(map);/*new double[s];
  if (map)
    memcpy(this->map, map, s*sizeof(double));
  else 
    memset(this->map, 0, s*sizeof(double));
  */
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
  path.pos = (int*)realloc(path.pos, count*3*sizeof(int));
  path.count = count;
  for (i = 0; i < count; ++i) {
    pos1 = (int*)cur->data;   
    memcpy(&path.pos[3*i], pos1, 3*sizeof(int));
    if (i > 0) {
      len += sqrt((pos1[0]-pos0[0])*(pos1[0]-pos0[0]) + (pos1[1]-pos0[1])*(pos1[1]-pos0[1]) + (pos1[2]-pos0[2])*(pos1[2]-pos0[2]));
    }
    pos0 = pos1;
    cur = cur->next;
  }
  path.len = len;
}

#define RAY_TRACE_STEP 0.1

void GridSearch3D::OptPath(const GridPath3D &path, GridPath3D &optPath) const
{
  double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, n, d;
  double dx0, dy0, dz0;
  double dx1, dy1, dz1;
  int i;
  int count = 1;
  double len = 0;
  int pos[3*path.count];

  if (path.len == 2) {
    optPath.pos = (int*)realloc(optPath.pos, 6*sizeof(int));
    memcpy(optPath.pos, path.pos, 6*sizeof(int));
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

  pos[0] = (int)x0;
  pos[1] = (int)y0;
  pos[2] = (int)z0;
  
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
	if (map[((int)(z0 + z))*length*width + ((int)(y0 + y))*length + (int)(x0 + x)]) {
	  pos[3*count] = (int)x1;
	  pos[3*count + 1] = (int)y1;
	  pos[3*count + 2] = (int)z1;
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
  
  pos[3*count] = (int)x2;
  pos[3*count + 1] = (int)y2;
  pos[3*count + 2] = (int)z2;
  len += sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));
  count++;
  optPath.pos = (int*)realloc(optPath.pos, count*3*sizeof(int));
  memcpy(optPath.pos, pos, count*3*sizeof(int));
  optPath.count = count;
  optPath.len = len;
}

