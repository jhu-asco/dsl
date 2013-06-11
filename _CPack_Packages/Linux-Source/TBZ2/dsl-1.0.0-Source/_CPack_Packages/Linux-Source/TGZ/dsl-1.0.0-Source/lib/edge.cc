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
