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
#include "travsearch.h"

using namespace std;
using namespace dsl;

TravSearch::TravSearch(int width, int height, const double *map, double scale) : 
  GridSearch(width, height, map, scale, false)
{
  // create all edges
  int i = 0;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x, ++i) {    
      for (int nbr = 0; nbr < NBR_COUNT; ++nbr) {
        int nx = x + NBR_OFFSETS[2*nbr];
        int ny = y + NBR_OFFSETS[2*nbr+1];
        if (nx < 0 || ny < 0 || nx >= width || ny >= height)
          continue;
        int ni = ny*width + nx;
        // create an edge b/n vertices at pos. i and ni
        Vertex *from = vertexMap[i];
        Vertex *to = vertexMap[ni];
      
        double ecost = CalcEdgeCost(this->map[i], this->map[ni], NBR_COSTS_[nbr]);
        //std::cout << ecost << std::endl;
        assert(ecost >= 0);
        Edge* edge = new Edge(from, to, ecost);
        
        graph.AddEdge(*edge);
      }
    }
  }
}

double TravSearch::CalcEdgeCost(double v1cost, double v2cost, double elength)
{
  assert(elength > 0);
  return scale*(fabs(v2cost- v1cost)/elength);
}

