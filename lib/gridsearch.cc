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
  


GridSearch::GridSearch(int width, int height, const GridCost& gridcost, const double *map, double scale) :
  Search(graph, gridcost),
  width(width),
  height(height), 
  scale(scale),
  cost(gridcost)
{
  int x, y, i, nbr, nx, ny, ni;

  int s = width*height;

  // initialize map
  this->map = new double[s];
  if (map)
    memcpy(this->map, map, s*sizeof(double));
  else 
    memset(this->map, 0, s*sizeof(double));

  // this will be used to map grid coordinates to dsl graph vertices
  vertexMap = new Cell2dVertex*[s];
  
  // create all vertices
  i = 0;
  for (y = 0; y < height; ++y) {
    for (x = 0; x < width; ++x, ++i) {
      //      VertexGridData* vgd = new VertexGridData();
      Cell2d cell(x, y, this->map[i]);
        //      vgd->p[0] = x;
        //      vgd->p[1] = y;
        //      vgd->cost = this->map[i];
      vertexMap[i] = new Vertex<Cell2d>(cell);
      //VertexGridData* vgdt = ( VertexGridData* )vertexMap[i]->data;
      //std::cout << vgdt->p[0] << " " << vgdt->p[1] << " " << vgdt->cost << std::endl;
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

        Cell2dVertex *from = vertexMap[i];
        Cell2dVertex *to = vertexMap[ni];
        
        // 
        //double ecost = edgeCost->CalcEdgeCost(this->map[i], this->map[ni], NBR_COSTS_[nbr]);
        //assert(ecost >= 0);
        double ecost = scale*gridcost.Real(from->data,to->data);
        Edge<Cell2d>* edge = new Edge<Cell2d>(from, to, ecost);
        //Edge* edge = new Edge(from, to, scale*ecost);
        

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
  Cell2dVertex *from = GetVertex(x1, y1);
  Cell2dVertex *to = GetVertex(x2, y2);
  if (!from || !to)
    return;
  
  double dx = x2-x1;
  double dy = y2-y1;
  Edge<Cell2d>* edge = new Edge<Cell2d>(from, to, scale*cost.Real(from->data, to->data)); 
  graph.AddEdge(*edge);
}


GridSearch::~GridSearch()
{
  for (int i = 0; i < width*height; ++i) {
    // delete[] (VertexGridData*)vertexMap[i]->data;
    //     vertexMap[i]->data = 0;
  }
  delete[] vertexMap;
  delete[] map;
}


double GridSearch::GetCost(int x, int y) const
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return 0;
  return map[y*width + x];
}


Cell2dVertex* GridSearch::GetVertex(int x, int y) const
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return 0;

  return vertexMap[y*width + x];
}

void GridSearch::RemoveVertex(int x, int y)
{
  Cell2dVertex *v = GetVertex(x, y);
  if (v)
    graph.RemoveVertex(*v);
}


void GridSearch::SetCost(int x, int y, double cost)
{
  std::map<int, Edge<Cell2d>*>::iterator ein, eout;

  
  int i = y*width + x;
 
  // if the cost has not changed simply return
  if (Eq(cost, this->map[i]))
    return;
  this->map[i] = cost;
  //  VertexGridData* v_data = (VertexGridData*)vertexMap[i]->data;
  Cell2d &v_data = vertexMap[i]->data;
  v_data.cost = cost;  

  // fix all connected edges
  ein = vertexMap[i]->in.begin();
  eout = vertexMap[i]->out.begin();
  for (;ein != vertexMap[i]->in.end(); ein++)
  {
    //int* pos = static_cast<int*>((*ein).second->from->data);

    //int dx = x-pos[0];
    //int dy = y-pos[1];
    
    //double ecost = edgeCost->CalcEdgeCost(GetCost(pos[0], pos[1]), cost, sqrt(dx*dx + dy*dy));
    //assert(ecost >= 0);
    ChangeCost(*ein->second, scale*this->cost.Real(ein->second->from->data, vertexMap[i]->data));
  }

  for (;eout != vertexMap[i]->out.end(); eout++)
  {
    //int* pos = static_cast<int*>((*eout).second->to->data);

    //int dx = x-pos[0];
    //int dy = y-pos[1];

    //double ecost = edgeCost->CalcEdgeCost(GetCost(pos[0], pos[1]), cost, sqrt(dx*dx + dy*dy));
    //assert(ecost >= 0);
    ChangeCost(*eout->second, scale*this->cost.Real(vertexMap[i]->data, (*eout).second->to->data));
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
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  Search::SetStart(*vertexMap[y*width + x]);
}

void GridSearch::SetGoal(int x, int y)
{
  if (x < 0 || x >= width || y < 0 || y >= height)
    return;
  Search::SetGoal(*vertexMap[y*width + x]);
}

void GridSearch::Plan(GridPath& path)
{
  double len = 0;
  Cell2dVertex *cur = start;
  int i;

  int count = Search<Cell2d>::Plan();
  path.cells.clear();
  //  path.pos = (int*)realloc(path.pos, count*2*sizeof(int));
  //  path.count = count;
  Cell2d prevCell;
  for (i = 0; i < count; ++i) {
    //    pos1 = ((VertexGridData*)cur->data)->p;   
    //Cell2d cell(cur->data);
    //    pos1 = cur->data.p;
    // memcpy(&path.pos[2*i], pos1, 2*sizeof(int));
    if (i > 0) {
      len += cur->data.Distance(prevCell);
    }
    prevCell = cur->data;
    path.cells.push_back(cur->data);
    cur = cur->next;
  }
  path.len = len;
}

#define RAY_TRACE_STEP 1.0

void GridSearch::OptPath(const GridPath &path, GridPath &optPath, double freeCost) const
{
  double x, y, x0, y0, x1, y1, x2, y2, n, d;
  double dx0, dy0;
  double dx1, dy1;
  int i;
  int count = 1;
  double len = 0;
  //  int pos[2*path.count];

  optPath.cells.clear();
  optPath.len = 0;

  if (path.len == 2) {
    optPath.cells = path.cells;
    optPath.len = path.len;
    return;
  }
  
  x2 = 0;
  y2 = 0;

  x0 = path.cells[0].p[0];
  y0 = path.cells[0].p[1];
  x1 = path.cells[1].p[0];
  y1 = path.cells[1].p[1];
  dx0 = x1 - x0;
  dy0 = y1 - y0;
  n = sqrt(dx0*dx0 + dy0*dy0);
  dx0 /= n;
  dy0 /= n;

  //  pos[0] = (int)x0;
  //  pos[1] = (int)y0;

  Cell2d cur((int)x0, (int)y0);
  optPath.cells.push_back(cur);
  
  for (i = 1; i < path.cells.size() - 1; ++i) {
    x1 = path.cells[i].p[0];
    y1 = path.cells[i].p[1];
    x2 = path.cells[i+1].p[0];
    y2 = path.cells[i+1].p[1];
    dx1 = x2 - x0;
    dy1 = y2 - y0;
    n = sqrt(dx1*dx1 + dy1*dy1);
    assert(n > .7);
    dx1 /= n;
    dy1 /= n;
    if (fabs((dx0-dx1)*(dx0-dx1) + (dy0-dy1)*(dy0-dy1)) > 1e-6) {
      n = sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2));
      for (d = RAY_TRACE_STEP; d < n; d += RAY_TRACE_STEP) {
	x = dx1*d;
	y = dy1*d;
	if (map[((int)(y0 + y))*width + (int)(x0 + x)] > freeCost) {
	  Cell2d cell((int)x1, (int)y1);
          optPath.cells.push_back(cell);
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
  
  Cell2d cell((int)x2, (int)y2);
  optPath.cells.push_back(cell);

  len += sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));
  optPath.len = len;
}
