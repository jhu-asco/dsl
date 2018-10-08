#include "carconnectivity.h"
#include <iostream>
#include "utils.h"


namespace dsl {

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using std::vector;


CarConnectivity::CarConnectivity(const CarGrid& grid,
                                 const vector<Eigen::Vector3d> &vs,
                                 double dt) : grid(grid), vs(vs), dt(dt)
{
}



CarConnectivity::CarConnectivity(const CarGrid& grid,
                                 double dt,
                                 double vx,
                                 double kmax,
                                 int kseg,
                                 bool onlyfwd)
    : grid(grid) {
  SetPrimitives(dt, vx, kmax, kseg, onlyfwd);
}

bool CarConnectivity::SetPrimitives(double dt, double vx, double kmax, int kseg, bool onlyfwd) {
  if (dt <= 0)
    return false;

  kseg = kseg < 1 ? 1 : kseg;

  this->dt = dt;

  vs.clear();

  for (int i = 0; i <= kseg; i++) {
    double k = i*kmax/kseg;
    double w = vx*k;
    vs.push_back(Vector3d(w, vx, 0));
    vs.push_back(Vector3d(-w, vx, 0));
    if (!onlyfwd) {
      vs.push_back(Vector3d(w, -vx, 0));
      vs.push_back(Vector3d(-w, -vx, 0));
    }
  }

  return true;
}

bool CarConnectivity::Flow(std::tuple< SE2Cell*, SE2Path, double>& pathTuple,
                           const Matrix3d& g0,
                           const Vector3d& v) const {

    //check if the cells encountered along the primitive and at the end, are free from obstacles or not
    double d = fabs(v[1]); // total distance along curve
    int n_seg = ceil(d/ (2 * grid.cs[1])); // 2 * grid.cs[1] is to improve efficiency
    double s = d/n_seg;
    SE2Cell* to = nullptr;
    SE2Path path;
    path.clear();
    for (int i_seg=1; i_seg<=n_seg; i_seg++) {
      Vector3d axy;
      Matrix3d g, dg;
      se2_exp(dg, (s*i_seg / d) * v);
      g = g0 * dg;
      se2_g2q(axy, g);
      to = grid.Get(axy);
      if (!to)
        return false;
      path.push_back(g);
    }

    pathTuple = std::make_tuple(to, path, d);
    return true;
}

bool CarConnectivity::
    operator()(const SE2Cell& from,
               std::vector< std::tuple<SE2Cell*, SE2Path, double> >& paths,
               bool fwd) const {
  Matrix3d g0;
  se2_q2g(g0, from.c);

  paths.clear();
  //  vector< Vector3d >::const_iterator it;
  for (const auto& s : vs) {
    // reverse time if fwd=false
    std::tuple<SE2Cell*, SE2Path, double> pathTuple;
    if (!Flow(pathTuple, g0, (fwd ? dt : -dt) * s))
      continue;

    assert(std::get<0>(pathTuple));


    // the path will now end inside the last cell but not exactly at the center,
    // which is a good enough
    // approximation if the cells are small

    // For exact trajectory generation to the center, we need to use more
    // complex inverse kinematics
    // which can be accomplished by uncommenting the following
    //    GenTraj(path.data, g0, path.cells.back().data, w,vx, vx, w,vx, dt/5);

    if (!std::get<1>(pathTuple).size())
      continue;

    // overwrite cost
    //    path.cost = cost; //.Real(path.cells.front(), path.cells.back());

    //    path.fwd = fwd;
    paths.push_back(pathTuple);
  }
  return true;
}


}
