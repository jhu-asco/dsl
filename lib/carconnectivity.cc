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

bool CarConnectivity::Flow(std::tuple< SE2CellPtr, SE2Path, double>& pathTuple,
                           const Matrix3d& g0,
                           const Vector3d& v, bool fwd) const {

    //check if the cells encountered along the primitive and at the end, are free from obstacles or not
    double d = fabs(v[1]); // total distance along curve
    int n_seg = ceil(d/ (2 * grid.cs[1])); // 2 * grid.cs[1] is to improve efficiency
    double s = d/n_seg;
    Vector3d axy;
    Matrix3d g, dg, g_inv;
    SE2CellPtr to(nullptr);
    for (int i_seg=1; i_seg<=n_seg; i_seg++) {
      se2_exp(dg, (s*i_seg / d) * v);
      g = g0 * dg;
      se2_g2q(axy, g);
      to = grid.Get(axy);
      if (!to)
        return false;
    }
    Matrix3d gf = to->data;

    //cost of the path
    if(fwd){
      se2_inv(g_inv,g0);
      dg = g_inv*gf;
    }else{
      se2_inv(g_inv,gf);
      dg = g_inv*g0;
    }
    Vector3d v_slip; se2_log(v_slip,dg);//twist that take you exactly to successor
    double cost = 0.1*abs(v_slip(0)) + v_slip.tail<2>().norm();//assuming carcost has ac=0.1

    //twist that takes you exactly only to (successor.x, successor.y) not angle
    Vector3d v_noslip;v_noslip << getWVx(dg(0,2),dg(1,2)),0;

    //If the graph expansion is in reverse, the waypoints are also added in reverse
    SE2Path path;
    path.clear();
    if(fwd){//when fwd g0 is the starting point
      for (int i_seg=1; i_seg<=n_seg; i_seg++) {
        Vector3d axy;
        Matrix3d g, dg;
        se2_exp(dg, (s*i_seg / d) * v_noslip);
        g = g0 * dg;
        se2_g2q(axy, g);
        path.push_back(g);
      }
    }else{//when fwd gf is the starting point
      for (int i_seg=1; i_seg<=n_seg; i_seg++) {
        Vector3d axy;
        Matrix3d g, dg;
        se2_exp(dg, (s*i_seg / d) * v_noslip);
        g = gf * dg;
        se2_g2q(axy, g);
        path.push_back(g);
      }
    }

    //Set the path tupule
    std::get<0>(pathTuple) = to;
    std::get<1>(pathTuple) = path;
    std::get<2>(pathTuple) = cost;

    return true;
}

bool CarConnectivity::
    operator()(const SE2Cell& from,
               std::vector< std::tuple<SE2CellPtr, SE2Path, double> >& paths,
               bool fwd) const {
  Matrix3d g0;
  se2_q2g(g0, from.c);

  paths.clear();
  //  vector< Vector3d >::const_iterator it;
  for (auto&& s : vs) {
    // reverse time if fwd=false
    std::tuple<SE2CellPtr, SE2Path, double> pathTuple;
    if (!Flow(pathTuple, g0, (fwd ? dt : -dt) * s,fwd))
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

bool CarConnectivity::GetPrims(const Vector3d pos, vector<vector<Vector2d>>& prims ){

  SE2CellPtr cell_start = grid.Get(pos);
  vector<std::tuple< SE2CellPtr, SE2Path, double> > paths;
  if(cell_start){
    (*this)(*cell_start,paths,true);
    Matrix3d g0 = cell_start->data;
    prims.reserve(paths.size());
    for(auto& path:paths){
      SE2CellPtr cell_to = std::get<0>(path);
      if(!cell_to)
        continue;
      SE2Path& gs=  std::get<1>(path);
      vector<Vector2d> prim(0); prim.reserve(gs.size()+1);
      prim.push_back(Vector2d(g0(0,2),g0(1,2)));
      for(auto& g: gs)
        prim.push_back(Vector2d(g(0,2),g(1,2)));
      prims.push_back(prim);
    }
    return true;
  }else{
    return false;
  }
}


/*
bool CarConnectivity::se2_times(Vector2d& t1,
               Vector2d& t2,
               Vector2d& t3,
               const Vector3d& s1,
               const Vector3d& s2,
               const Vector3d& s3,
               const Matrix3d& g) {
  const double& w1 = s1(0);
  const double& v1 = s1(1);
  const double& w2 = s2(0);
  const double& v2 = s2(1);
  const double& w3 = s3(0);
  const double& v3 = s3(1);
  const double& xf = g(0, 2);
  const double& yf = g(1, 2);
  const double& ca = g(0, 0);
  const double& sa = g(1, 0);

  double s = (w1 * w3 * (2 * v1 * v3 - 2 * v1 * v3 * ca + w1 * w3 * xf * xf +
                         w1 * w3 * yf * yf - 2 * v1 * w3 * yf +
                         2 * v3 * w1 * yf * ca - 2 * v3 * w1 * xf * sa));
  if (s < 0)
    return false;

  double ss = sqrt(s);

  t2[0] = ss / (v2 * w1 * w3);
  t2[1] = -t2[0];

  double d = (v3 * w1 - 2 * v1 * w3 + v3 * w1 * ca + w1 * w3 * yf);
  if (d < 1e-32)
    return false;

  double n1 = v3 * w1 * sa - w1 * w3 * xf;
  double k1[2] = {(ss + n1) / d, (-ss + n1) / d};

  double af = atan2(sa, ca);

  t1[0] = 2 * atan(k1[0]) / w1;
  t1[1] = 2 * atan(k1[1]) / w1;

  t3[0] = fangle(af - t1[0] * w1) / w3;
  t3[1] = fangle(af - t1[1] * w1) / w3;

  return true;
}

static double CarConnectivity::GenTraj(vector< Matrix3d >& gs,
                      const Matrix3d& g0,
                      const Vector3d& s,
                      double t0,
                      double tf,
                      double dt) {
  if (t0 >= tf)
    return t0;

  Matrix3d g;
  double t = t0;
  for (; t <= tf; t += dt) {
    se2_exp(g, t * s);
    gs.push_back(g0 * g);
  }
  return t;
}

bool CarConnectivity::GenTraj(vector< Matrix3d > gs,
             const Matrix3d& g0,
             const Matrix3d& gf,
             double w1,
             double v1,
             double v2,
             double w3,
             double v3,
             double dt,
             double eps = 1e-16) {
  Matrix3d gi;
  se2_inv(gi, g0);
  Matrix3d dg = gi * gf; // relative pose

  Vector2d ts[3]; // each of the three segments has two solutions

  // four possible sequences of three vector fields
  Vector3d ss[4][3] = {
      {Vector3d(w1, v1, 0), Vector3d(0, v2, 0), Vector3d(w3, v3, 0)},
      {Vector3d(-w1, v1, 0), Vector3d(0, v2, 0), Vector3d(w3, v3, 0)},
      {Vector3d(w1, v1, 0), Vector3d(0, v2, 0), Vector3d(-w3, v3, 0)},
      {Vector3d(-w1, v1, 0), Vector3d(0, v2, 0), Vector3d(-w3, v3, 0)}};

  double tmin = 1e10; // best total time
  int ib = -1;        // best index i
  int ij = -1;        // best index j
  double tbs[3];      // best times
  for (int i = 0; i < 4; ++i) {
    bool ok = se2_times(ts[0], ts[1], ts[2], ss[i][0], ss[i][1], ss[i][2], dg);
    if (!ok)
      continue;
    for (int j = 0; j < 2; ++j) { // there are two time choices
      double ttot = 0;
      for (int k = 0; k < 3; ++k)
        ttot += fabs(ts[k][j]); // compute total time for this choice

      if (ttot < tmin) {
        tmin = ttot;
        ib = i;
        ij = j;
        tbs[0] = ts[0][j];
        tbs[1] = ts[1][j];
        tbs[2] = ts[2][j];
      }
    }
  }
  assert(ib >= 0);
  assert(ij >= 0);

  // integrate the best sequence forward
  Matrix3d gr = g0; // reference
  double t0 = 0;
  for (int k = 0; k < 3; ++k) {
    Vector3d& s = ss[ib][k];
    double& tb = tbs[k];
    if (fabs(tb) < eps)
      continue;

    if (tb < 0) {
      tb = -tb;
      s = -s;
    }

    //    cout << "k=" << k << " tb=" << tb << " s=" << s.transpose() << endl;

    vector< Matrix3d > gis;

    // if the next segment should contains a sample
    //    if (tb - t0 >= dt) {
    double tl = GenTraj(gis,
                        gr,
                        s,
                        t0,
                        tb,
                        dt); // returns time of last pose, tl<=tb & (tb-tl)<dt

    gs.insert(gs.end(), gis.begin(), gis.end()); // add to path
    t0 = dt + tl -
        tb; // update start time so that we have even time-steps across sections
    //    } else {
    //      t0 =
    //    }

    se2_exp(dg, tb * s); // total displacement along this segment
    gr = gr * dg;        // update start of next section

    //    cout << "dt=" << dt << " tl=" << tl << " t0=" << t0 << endl;
  }

  return true;
}
*/

}
