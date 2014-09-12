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
#include <set>
#include <string.h>
#include "gridsearch3dvel.h"

using namespace std;
using namespace dsl;

#define _USE_MATH_DEFINES
#define MAX(a,b) (a>b?a:b)


GridSearch3DVel::GridSearch3DVel(int length, int width, int height, int numYaws, int numPitches, const double *map, double scale) : 
  Search(graph, cost),
  length(length),
  width(width),
  height(height), 
  numYaws(numYaws),
  numPitches(numPitches),
  scale(scale)
{
  int s = length*width*height;
  this->map = new double[s];
  if (map)
    memcpy(this->map, map, s*sizeof(double));
  else 
    memset(this->map, 0, s*sizeof(double));
}

void GridSearch3DVel::Init()
{  
  int x, y, z, theta, phi, i, nx, ny, nz;
  int* pos;

  int s = length*width*height;
  double degPerYaw = 360./numYaws;
  double degPerPitch = 360./numPitches;

  // initialize map
  //this->map = const_cast<double*>(map);/*
  
  cout << "Grid map initialized" << endl;
  // this will be used to map grid coordinates to dsl graph vertices
  vertexMap = new Vertex*[s*numYaws*numPitches];
  cout << "Num Vertices: " << s*numYaws*numPitches << endl;  
  // create all vertices
  i = 0;
  for(z = 0; z < height; ++z){
    for (y = 0; y < width; ++y) {
      for (x = 0; x < length; ++x) {
        for(theta = 0; theta < numYaws; theta++){
          for(phi = 0; phi < numPitches; phi++, i++){
            pos = new int[5];
            pos[0] = x;
            pos[1] = y;
            pos[2] = z;
            pos[3] = (int)theta*degPerYaw; 
            pos[4] = (int)phi*degPerPitch + 90; 
            vertexMap[i] = new Vertex(pos);
            graph.AddVertex(*vertexMap[i]);
          }
        }
      }     
    }
  }

  cout << "Created vertices" << endl;
  // create all edges by trying to connect each state to a state that lies on some sphere around it
  i = 0;
  int lastperc = 0;
  int costi = 0;
  for(z = 0; z < height; ++z){
    for (y = 0; y < width; ++y) {
      for (x = 0; x < length; ++x, costi++) {    
        for(theta = 0; theta < numYaws; theta++){
          for(phi = 0; phi < numPitches; phi++, i++){
            int perc = ((100*i)/(s*numYaws*numPitches));
            if(perc != lastperc)
            {
              cout << perc << endl;
              lastperc = perc;
            }
            if((int)(this->map[z*width*length + y*length + x]) == DSL3D_OCCUPIED)
              continue;
	    Vertex *from = vertexMap[i];
	    //int* pos = (int*)from->data;
            //cout << pos[0] << " " << pos[1] << " " << pos[2] << endl;
            
            std::set<int> connectedCells;    
            // connect to sphere around point
            for(int st = 0; st < 360; st+=30){
              for(int sp = 0; sp < 360; sp+=30){
                for(int r = 1; r < 3; r++){
                  double trad = st*M_PI/180.;
                  double prad = sp*M_PI/180.;
                  nx = x + r*cos(trad)*sin(prad);  
                  ny = y + r*sin(trad)*sin(prad); 
                  nz = z + r*cos(prad);
                  if(nx < 0 || ny < 0 || nz < 0 || nx >= length || ny >= width || nz >= height)
                    continue;
                 
                  int ncosti = nz*width*length + ny*length + nx;

                  // Don't connect to same cell twice
                  if(connectedCells.find(ncosti) != connectedCells.end())
                    continue;
                  connectedCells.insert(ncosti);

                  // Don't connect to occupied cell
                  if((int)(this->map[ncosti]) == DSL3D_OCCUPIED)
                    continue;

                  for(int ht = 0; ht < numYaws; ht++){
                    for(int hp = 0; hp < numPitches; hp++){
                      Vertex *to = GetVertex(nx,ny,nz,ht,hp);
	              //int* pos = (int*)to->data;
                      //cout << "to: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
                      
                      double cost = -1;
                      if((int)(cost = GetPathCost(*from, *to)) != -1)      
                      {  
                        Edge* edge = new Edge(from, to, scale*(MAX(this->map[costi], this->map[ncosti]) + cost));
                        //cout << "cost: " << (scale*(MAX(this->map[costi], this->map[ncosti]) + cost)) << endl;
                        graph.AddEdge(*edge);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    //cout << z << endl;
  }
}

double GridSearch3DVel::GetPathCost(const Vertex &from, const Vertex &to) const
{
  GridPath3D path;
  GetTrajectory(from, to, path);
  return path.len;
}

void GridSearch3DVel::AddEdge(int x1, int y1, int z1, int t1, int p1, int x2, int y2, int z2, int t2, int p2)
{  
  Vertex *from = GetVertex(x1, y1, z1, t1, p1);
  Vertex *to = GetVertex(x2, y2, z2, t2, p2);
  if (!from || !to)
    return;
  
  Edge* edge = new Edge(from, to, scale*(MAX(GetCost(x1, y1, z1), GetCost(x2, y2, z2)) + GetPathCost(from, to))); 
  graph.AddEdge(*edge);
}


GridSearch3DVel::~GridSearch3DVel()
{
  for (int i = 0; i < length*width*height*numYaws*numPitches; ++i) {
    delete[] (int*)vertexMap[i]->data;
    vertexMap[i]->data = 0;
  }
  delete[] vertexMap;
  delete[] map;
}


double GridSearch3DVel::GetCost(int x, int y, int z) const
{
  return map[z*length*width + y*length + x];
}


Vertex* GridSearch3DVel::GetVertex(int x, int y, int z, int t, int p) const
{
  if (x < 0 || x >= length || y < 0 || y >= width || z < 0 || z >= height || t < 0 || t >= numYaws || p < 0 || p >= numPitches)
    return 0;

  return vertexMap[z*length*width*numYaws*numPitches + y*length*numYaws*numPitches + x*numYaws*numPitches + t*numPitches + p];
}

void GridSearch3DVel::RemoveVertex(int x, int y, int z, int t, int p)
{
  Vertex *v = GetVertex(x, y, z, t, p);
  if (v)
    graph.RemoveVertex(*v);
}


void GridSearch3DVel::SetCost(int x, int y, int z, int t, int p, double cost)
{
  std::map<int, Edge*>::iterator ein, eout;

  assert(cost >= 0);
  
  int mapi = z*length*width + y*length + x;
  // if the cost has not changed simply return
  if (Eq(cost, this->map[mapi]))
    return;
  this->map[mapi] = cost;
  
  // fix all connected edges
  Vertex* v = GetVertex(x,y,z,t,p); 
  ein = v->in.begin();
  eout = v->out.begin();
  for (;ein != v->in.end(); ein++)
    ChangeCost(*ein->second, cost > 10 ? 10000 : 0);

  for (;eout != v->out.end(); eout++)
    ChangeCost(*eout->second, cost > 10 ? 10000 : 0);
}

void GridSearch3DVel::SetStart(int x, int y, int z, int t, int p)
{
  Search::SetStart(*(GetVertex(x,y,z,t,p)));
}

void GridSearch3DVel::SetGoal(int x, int y, int z, int t, int p)
{
  Search::SetGoal(*(GetVertex(x,y,z,t,p)));
}

void GridSearch3DVel::Plan(GridPath3D& path)
{
  path = GridPath3D();
  Vertex *last = start, *cur;
  int i;
  int count;

  count = Search::Plan();
  cur = start->next;
  cout << count << endl;
  //path.pos = (double*)realloc(path.pos, count*3*sizeof(double));
  //path.count = count;
  for (i = 1; i < count; ++i) {
    GridPath3D nextPath;
    GetTrajectory(*last,*cur, nextPath);
    path.AppendPath(nextPath);
    last = cur;
    cur = cur->next;
  }
}
