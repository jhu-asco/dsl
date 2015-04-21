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
#include "gridsearch.h"

using namespace std;
using namespace dsl;

#define MAX(a,b) (a>b?a:b)
// cell neighbor connectivity
// this can be either 4 (faster,less accurate) or 8 (slower, more accurate) 
#define NBR_COUNT 8
#define SQRT2 1.414213562373095

const int GridSearch::NBR_OFFSETS[8*2] = {-1,-1, 0,-1, 1,-1, -1,0, 1,0, -1,1, 0,1 ,1,1};
const double GridSearch::NBR_COSTS_[8] = {SQRT2,  1.0, SQRT2, 1.0, 1.0, SQRT2, 1.0, SQRT2};

GridSearch::GridSearch(int width, int height, EdgeCost* edgeCost, const double *map, double scale) : 
  Search(graph, cost),
  width(width),
  height(height), 
  scale(scale),
  edgeCost(edgeCost)
{
  int x, y, i, nbr, nx, ny, ni;
  int* pos;

  int s = width*height;

  // initialize map
  this->map = new double[s];
  if (map)
    memcpy(this->map, map, s*sizeof(double));
  else 
    memset(this->map, 0, s*sizeof(double));

  // this will be used to map grid coordinates to dsl graph vertices
  vertexMap = new Vertex*[s];
  
  // create all vertices
  i = 0;
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x, ++i) {
      pos = new int[2];
      pos[0] = x;
      pos[1] = y;
      vertexMap[i] = new Vertex(pos);
      graph.AddVertex(*vertexMap[i]);
                  //      this->vertices.push_back();
    }
  }


  // create all edges
  i = 0;
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x, ++i) {    
      for (nbr = 0; nbr < NBR_COUNT; ++nbr) {
        nx = x + NBR_OFFSETS[2*nbr];
        ny = y + NBR_OFFSETS[2*nbr+1];
        if (nx < 0 || ny < 0 || nx >= width || ny >= height)
          continue;
        ni = ny*width + nx;
        // create an edge b/n vertices at pos. i and ni

        Vertex *from = vertexMap[i];
        Vertex *to = vertexMap[ni];
        
        double ecost = edgeCost->CalcEdgeCost(this->map[i], this->map[ni], NBR_COSTS_[nbr]);
        assert(ecost >= 0);
        //Edge* edge = new Edge(from, to, scale*(MAX(this->map[i], this->map[ni]) + NBR_COSTS_[nbr]));
        Edge* edge = new Edge(from, to, scale*ecost);
        

        graph.AddEdge(*edge);
        //        from->out[edge->id] = edge;
        //        to->in[edge->id] = edge;
        //        edges.push_back(edge);
        
      }
    }
  }
}

void GridSearch::AddEdge(int x1, int y1, int x2, int y2)
{  
  Vertex *from = GetVertex(x1, y1);
  Vertex *to = GetVertex(x2, y2);
  if (!from || !to)
    return;
  
  double dx = x2-x1;
  double dy = y2-y1;
  Edge* edge = new Edge(from, to, sqrt(dx*dx + dy*dy)); 
  graph.AddEdge(*edge);
}


GridSearch::~GridSearch()
{
  for (int i = 0; i < width*height; ++i) {
    delete[] (int*)vertexMap[i]->data;
    vertexMap[i]->data = 0;
  }
  delete[] vertexMap;
  delete[] map;
}


double GridSearch::GetCost(int x, int y) const
{
  return map[y*width + x];
}


Vertex* GridSearch::GetVertex(int x, int y) const
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return 0;

  return vertexMap[y*width + x];
}

void GridSearch::RemoveVertex(int x, int y)
{
  Vertex *v = GetVertex(x, y);
  if (v)
    graph.RemoveVertex(*v);
}


void GridSearch::SetCost(int x, int y, double cost)
{
  std::map<int, Edge*>::iterator ein, eout;

  
  int i = y*width + x;
 
  // if the cost has not changed simply return
  if (Eq(cost, this->map[i]))
    return;
  this->map[i] = cost;
  
  // fix all connected edges
  ein = vertexMap[i]->in.begin();
  eout = vertexMap[i]->out.begin();
  for (;ein != vertexMap[i]->in.end(); ein++)
  {
    int* pos = static_cast<int*>((*ein).second->from->data);

    int dx = x-pos[0];
    int dy = y-pos[1];
    
    double ecost = edgeCost->CalcEdgeCost(GetCost(pos[0], pos[1]), cost, sqrt(dx*dx + dy*dy));
    assert(ecost >= 0);
    ChangeCost(*ein->second, scale*ecost);
  }

  for (;eout != vertexMap[i]->out.end(); eout++)
  {
    int* pos = static_cast<int*>((*eout).second->to->data);

    int dx = x-pos[0];
    int dy = y-pos[1];

    double ecost = edgeCost->CalcEdgeCost(GetCost(pos[0], pos[1]), cost, sqrt(dx*dx + dy*dy));
    assert(ecost >= 0);
    ChangeCost(*eout->second, scale*ecost);
  }
  //for (;ein != vertexMap[i]->in.end(); ein++)
  //  ChangeCost(*ein->second, cost > 10 ? 10000 : 0);

  //for (;eout != vertexMap[i]->out.end(); eout++)
  //  ChangeCost(*eout->second, cost > 10 ? 10000 : 0);
}

void GridSearch::SetMap(const double *map)
{
  int x, y, i;
  i = 0;
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x, ++i) {
      SetCost(x, y, map[i]);
    }
  }
  memcpy(this->map, map, width*height*sizeof(double));
}


void GridSearch::SetStart(int x, int y)
{
  Search::SetStart(*vertexMap[y*width + x]);
}

void GridSearch::SetGoal(int x, int y)
{
  Search::SetGoal(*vertexMap[y*width + x]);
}

void GridSearch::Plan(GridPath& path)
{
  int count;
  double len = 0;
  Vertex *cur = start;
  int *pos1 = 0, *pos0 = 0;
  int i;

  count = Search::Plan();
  path.pos = (int*)realloc(path.pos, count*2*sizeof(int));
  path.count = count;
  for (i = 0; i < count; ++i) {
    pos1 = (int*)cur->data;   
    memcpy(&path.pos[2*i], pos1, 2*sizeof(int));
    if (i > 0) {
      len += sqrt((pos1[0]-pos0[0])*(pos1[0]-pos0[0]) + (pos1[1]-pos0[1])*(pos1[1]-pos0[1]));
    }
    pos0 = pos1;
    cur = cur->next;
  }
  path.len = len;
}

#define RAY_TRACE_STEP 1.0

void GridSearch::OptPath(const GridPath &path, GridPath &optPath) const
{
  double x, y, x0, y0, x1, y1, x2, y2, n, d;
  double dx0, dy0;
  double dx1, dy1;
  int i;
  int count = 1;
  double len = 0;
  int pos[2*path.count];

  if (path.len == 2) {
    optPath.pos = (int*)realloc(optPath.pos, 4*sizeof(int));
    memcpy(optPath.pos, path.pos, 4*sizeof(int));
    optPath.count = 2;
    optPath.len = path.len;
    return;
  }
  
  x2 = 0;
  y2 = 0;

  x0 = path.pos[0];
  y0 = path.pos[1];
  x1 = path.pos[2];
  y1 = path.pos[3];
  dx0 = x1 - x0;
  dy0 = y1 - y0;
  n = sqrt(dx0*dx0 + dy0*dy0);
  dx0 /= n;
  dy0 /= n;

  pos[0] = (int)x0;
  pos[1] = (int)y0;
  
  for (i = 1; i < path.count - 1; ++i) {
    x1 = path.pos[2*i];
    y1 = path.pos[2*i + 1];
    x2 = path.pos[2*(i+1)];
    y2 = path.pos[2*(i+1)+1];
    dx1 = x2 - x0;
    dy1 = y2 - y0;
    n = sqrt(dx1*dx1 + dy1*dy1);
    assert(n > .7);
    dx1 /= n;
    dy1 /= n;
    if (fabs((dx0-dx1)*(dx0-dx1) + (dy0-dy1)*(dy0-dy1)) > .00001) {
      n = sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2));
      for (d = RAY_TRACE_STEP; d < n; d += RAY_TRACE_STEP) {
	x = dx1*d;
	y = dy1*d;
	if (map[((int)(y0 + y))*width + (int)(x0 + x)]) {
	  pos[2*count] = (int)x1;
	  pos[2*count + 1] = (int)y1;
	  count++;
	  len += sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
	  x0 = x1;
	  y0 = y1;
	  break;
	}
      }
      dx0 = x2 - x0;
      dy0 = y2 - y0;
      n = sqrt(dx0*dx0 + dy0*dy0);
      dx0 /= n;
      dy0 /= n;
    }
  }
  
  pos[2*count] = (int)x2;
  pos[2*count + 1] = (int)y2;
  len += sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));
  count++;
  optPath.pos = (int*)realloc(optPath.pos, count*2*sizeof(int));
  memcpy(optPath.pos, pos, count*2*sizeof(int));
  optPath.count = count;
  optPath.len = len;
}

