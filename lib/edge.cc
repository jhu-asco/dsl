// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "edge.h"
#include "vertex.h"

using namespace dsl;

int Edge::s_id = 0;


Edge::Edge(Vertex *from, Vertex *to, double cost) :
  id(s_id), from(from), to(to), cost(cost), costChange(-INF) 
{
  ++s_id;
}    

Edge::~Edge()
{
}

namespace dsl {
  std::ostream& operator<<(std::ostream &os, const Edge &e)
  {
    os << e.id << " ";
    if (e.from)
      os << e.from->id << " ";
    else 
      os << "-1 ";
    if (e.to)
      os << e.to->id << " ";
    else 
      os << "-1 ";
    os << e.cost << " ";
    os << e.costChange;
    return os;
  };
  
  std::istream& operator>>(std::istream &is, Edge &e)
  {   
    return is;
  }
}
