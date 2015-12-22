#include "carcost.h"

using namespace dsl;

static double fangle(double a) {
  double d = fmod(a, 2*M_PI);
  if (d > M_PI)
    d -= 2*M_PI;
  else
    if (d < -M_PI)
      d += 2*M_PI;
  // theta should now be in [-pi,pi]
  return d;
}


double CarCost::Real(const Cell<3> &a, const Cell<3> &b) const {
  // default real cost is euclidean distance + average cell cost multiplied by Euclidean distance
  double d1 = (a.c.head<2>() - b.c.head<2>()).norm(); // position distance
  double d2 = fabs(fangle(a.c(2) - b.c(2)));          // orientation distance
  return (1 + (a.cost + b.cost)/2)*(d1 + d2);
}

double CarCost::Heur(const Cell<3> &a, const Cell<3> &b) const {
  // default Heuristic cost is the Euclidean distance

  double d1 = (a.c.head<2>() - b.c.head<2>()).norm(); // position distance
  double d2 = fabs(fangle(a.c(2) - b.c(2)));          // orientation distance
  
  return d1 + d2;
}

