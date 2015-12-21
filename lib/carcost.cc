#include "carcost.h"

using namespace dsl;

double CarCost::Real(const Cell<3> &a, const Cell<3> &b) const {
  // default real cost is euclidean distance + average cell cost multiplied by Euclidean distance
  return (1 + (a.cost + b.cost)/2)*(a.c.head<2>() - b.c.head<2>()).norm();
}

double CarCost::Heur(const Cell<3> &a, const Cell<3> &b) const {
  // default Heuristic cost is the Euclidean distance
  return (a.c.head<2>() - b.c.head<2>()).norm();
}
