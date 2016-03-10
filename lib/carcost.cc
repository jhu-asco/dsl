#include "carcost.h"
#include "utils.h"

namespace dsl {

CarCost::CarCost(double ac, double eps) : ac(ac), eps(eps) {
  assert(ac >= 0);
  assert(eps > 0 && eps < 1);
}

double CarCost::Real(const SE2Cell& a, const SE2Cell& b) const {
  double d1 = (a.c.tail< 2 >() - b.c.tail< 2 >()).norm(); // position distance
  double d2 = fabs(fangle(a.c(0) - b.c(0))); // orientation distance
  return d1 + ac * d2;
}

double CarCost::Heur(const SE2Cell& a, const SE2Cell& b) const {
  double d1 = (a.c.tail< 2 >() - b.c.tail< 2 >()).norm(); // position distance
  double d2 = fabs(fangle(a.c[0] - b.c[0])); // orientation distance

  return (1 - eps)*(d1 + ac * d2);
}

}
