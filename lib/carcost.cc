#include "carcost.h"
#include "utils.h"

using namespace dsl;

CarCost::CarCost(double ac) : ac(ac)
{
}

double CarCost::Real(const SE2Cell &a, const SE2Cell &b) const {
  // default real cost is euclidean distance + average cell cost multiplied by Euclidean distance
  double d1 = (a.c.tail<2>() - b.c.tail<2>()).norm(); // position distance
  double d2 = fabs(fangle(a.c(0) - b.c(0)));          // orientation distance
  return (1 + (a.cost + b.cost)/2)*(d1 + ac*d2);
}

double CarCost::Heur(const SE2Cell &a, const SE2Cell &b) const {
  // default Heuristic cost is the Euclidean distance
  
  return Real(a,b);
  //  double d1 = (a.c.tail<2>() - b.c.tail<2>()).norm(); // position distance
  //  double d2 = fabs(fangle(a.c(0) - b.c(0)));          // orientation distance
  
  //  return d1 + ac*d2;
}

