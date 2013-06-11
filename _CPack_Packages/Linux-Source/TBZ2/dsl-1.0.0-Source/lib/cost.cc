#include "cost.h"

using namespace dsl;

double Cost::Real(const Vertex &va, const Vertex &vb) const
{
  return Heur(va, vb) + 1e-10;
}
