#include "gridcost.h"
#include <cmath>

using namespace dsl;


double GridCost::Real(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = s1pos[0] - s2pos[0];
  double dy = s1pos[1] - s2pos[1];
  return sqrt(dx*dx + dy*dy);
}


double GridCost::Heur(const Vertex &a, const Vertex &b) const
{
  int* s1pos = (int*)a.data;
  int* s2pos = (int*)b.data;
  double dx = fabs((double)(s1pos[0] - s2pos[0]));
  double dy = fabs((double)(s1pos[1] - s2pos[1]));
  if (dx > dy)
    return dx;
  else
    return dy;
}
