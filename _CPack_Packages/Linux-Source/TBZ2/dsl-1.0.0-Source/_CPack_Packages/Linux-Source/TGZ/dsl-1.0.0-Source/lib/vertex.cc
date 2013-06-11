#include "vertex.h"
#include "edge.h"


using namespace dsl;

int Vertex::s_id = 0;

Vertex::Vertex(void *data) : 
  id(s_id), data(data)
{
  Reset();
  ++s_id;
}

Vertex::~Vertex()
{
  in.clear();
  out.clear();
}

void Vertex::Reset()
{
  next = 0;
  prev = 0;
  rhs = g = INF;
  t = Vertex::NEW;
  openListNode = 0;
  key[0] = key[1] = INF;
  
#ifdef DSL_DSTAR_MIN
  r = 0;
#endif
}

Edge* Vertex::Find(const Vertex &v, bool in) const
{
  std::map<int, Edge*>::const_iterator it;
  
  if (in)
    for (it = this->in.begin(); it != this->in.end(); ++it) {
      Edge *e = it->second;
      if (e->from == &v)
        return e;
    }
  else
    for (it = this->out.begin(); it != this->out.end(); ++it) {
      Edge *e = it->second;
      if (e->to == &v)
        return e;
    }
  return 0;
}

namespace dsl {
  std::ostream& operator<<(std::ostream &os, const Vertex &v)
  {
    os << v.id << " ";
    std::map<int, Edge*>::const_iterator it;
    os << "[";
    for (it = v.in.begin(); it != v.in.end(); ++it) {
      os << it->second->id << " ";
    }
    os << "]<-*->[";
    for (it = v.out.begin(); it != v.out.end(); ++it) {
      os << it->second->id << " ";
    }
    os << "] ";

    if (v.next)
      os << v.next->id << " ";
    else 
      os << "-1 ";
    
    os << "rhs=" << v.rhs << " g=" << v.g << " t=" << v.t << " node=" << v.openListNode << " key=(" << v.key[0] << "," << v.key[1] << ")" << std::endl;;
    

    return os;
  };
  
  std::istream& operator>>(std::istream &is, Vertex &v)
  {   
    return is;
  }
}
