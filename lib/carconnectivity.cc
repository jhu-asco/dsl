#include "carconnectivity.h"
#include <iostream>

using namespace dsl;
using namespace std;


static double DubinsTimes(double ts[3], 
                          double ws[2],
                          double v,
                          double w,
                          Matrix3d g,
                          double tmin = 1e-6)
{
  // tb = -((v*sin(ta*wa))/wa - x + (v*(st - sin(ta*wa)))/wb)/(v*cos(ta*wa))
  //  ((2*v*wb - wa*v - ct*wa*v - wa*wb*y)*k^2 + (2*st*wa*v - 2*wa*wb*x)*k + ct*wa*v - wa*v + wa*wb*y)/((wa*wb)*k^2 - wa*wb)
  
  const double& ct = g(0,0);
  const double& st = g(1,0);
  const double& x = g(0,2);
  const double& y = g(1,2);
  
  double was[4] = {w, w, -w, -w};
  double wcs[4] = {w, -w, w, -w};
  
  double ttb = INFINITY;  // total time best
  bool ok = false;
  
  for (int j = 0; j < 4; ++j) {  
    
    double &wa = was[j];
    double &wc = wcs[j];

    double a = (2*v*wc - wa*v - ct*wa*v - wa*wc*y);
    double b = (2*st*wa*v - 2*wa*wc*x);
    double c = ct*wa*v - wa*v + wa*wc*y;
    
    double s = b*b - 4*a*c;
    //  cout << "s=" << s << " a=" << a << endl;
    
    if (s < 0)
      return -1;
    
    s = sqrt(s);
    
    double k[2] = { (-b + s)/(2*a), (-b - s)/(2*a)};
        
    // go through all solutions
    for (int i = 0; i < 2; ++i) {
      
      double k2 = k[i]*k[i];
      double cta = (1 - k2)/(1 + k2);
      double sta = 2*k[i]/(1 + k2);
      
      double ctc = ct*cta + st*sta;
      double stc = st*cta - ct*sta;
      
      double tbt = -((v*sta)/wa - x + (v*(st - sta))/wc)/(v*cta);
      
      if (tbt < tmin)
        continue;
      
      double tct = atan2(stc, ctc)/wc;
      if (tct < 0)
        tct += ((wc < 0 ? -2*M_PI : 2*M_PI)/wc);
      assert(tct > -tmin);
      
      double tat = atan2(sta, cta)/wa;
      if (tat < 0)
        tat += ((wa < 0 ? -2*M_PI : 2*M_PI)/wa);
      assert(tat > -tmin);
      
      double tt = tat + tbt + tct;  // total time
      if (tt < ttb) {
        ts[0] = tat;
        ts[1] = tbt;
        ts[2] = tct;
        ws[0] = wa;
        ws[1] = wc;
        ttb = tt;
        ok = true;
      }
    }
  }

  if (ok) {
    return ttb;
  }
    
  return -1;
}


static void q2g(Matrix3d &m, const Vector3d &q)
{
  double ct = cos(q[2]);
  double st = sin(q[2]);

  m(0,0) = ct;   m(0,1) = -st;  m(0,2) = q[0];
  m(1,0) = st;   m(1,1) = ct;   m(1,2) = q[1];
  m(2,0) = 0;    m(2,1) = 0;    m(2,2) = 1;
}

static void g2q(Vector3d &q, const Matrix3d &m) 
{
  q[0] = m(0,2);
  q[1] = m(1,2);
  q[2] = atan2(m(1,0), m(0,0));
}

static void exp(Matrix3d &m, const Vector3d &v, double tol = 1e-16)
{
  const double &w = v[2];

  if (fabs(w) < tol) {
    m(0,0) = 1; m(0,1) = 0; m(0,2) = v(0);
    m(1,0) = 0; m(1,1) = 1; m(1,2) = v(1);
    m(2,0) = 0; m(2,1) = 0; m(2,2) = 1;
    return;
  }
  
  double c = cos(w);
  double s = sin(w);
  
  double ax = v[1]/w;
  double ay = -v[0]/w;

  m(0,0) = c; m(0,1) = -s; m(0,2) = (c - 1)*ax - s*ay;
  m(1,0) = s; m(1,1) = c;  m(1,2) = s*ax + (c - 1)*ay;
  m(2,0) = 0; m(2,1) = 0;  m(2,2) = 1;
}


CarConnectivity::CarConnectivity(const CarGrid &grid) : grid(grid) 
{
  double tphi = 0.577; // tan(M_PI/6); // steering angle
  v = 1;
  SetPrimitives(v, tphi*v, 1);
}

bool CarConnectivity::SetPrimitives(double v, double w, double dt) {
  if (dt <= 0)
    return false;
    
  this->v = v;
  this->w = w;
  this->dt = dt;
  vs.push_back(Vector3d(v, 0, 0));
  vs.push_back(Vector3d(v, 0, w));
  vs.push_back(Vector3d(v, 0, -w));
  vs.push_back(Vector3d(-v, 0, 0));
  vs.push_back(Vector3d(-v, 0, w));
  vs.push_back(Vector3d(-v, 0, -w));

  return false;
}


bool CarConnectivity::Flow(GridPath<3>& path, const Matrix3d &g0, const Vector3d &v) const
{
  double d = fabs(v[0]);
  double s = 2*grid.cs[0];  // set step-size to side-length
  
  Matrix3d g;
  Vector3d q;

  //  cout << "d=" << d << " s=" << s << endl;

  path.cells.clear();
  path.len = 0;

  // generate it backwards to more efficiently handle obstacles
  //  for (double a = d; a > 0; a -= s) {
  for (double a = s; a <= d; a += s) {
    Matrix3d dg;
    exp(dg, (a/d)*v);
    g2q(q, g0*dg);
    
    if (!grid.Valid(q)) {
      return false;
    }
    
    int id = grid.Id(q);
    Cell<3> *cell = grid.cells[id];
    if (!cell || cell->cost > grid.maxCost) {
      return false;      
    }

    path.cells.push_back(*cell);
    path.len += cell->cost;   // add up all cost along cells
  }
  path.len = (1 + path.len)*fabs(v[0]);   // regard cell cost as "traversability" which additionally penalizes the travelled distance

  return true;
}




bool CarConnectivity::operator()(const Cell<3>& from, 
                                 vector<GridPath<3> >& paths, 
                                 bool fwd) const {

  Matrix3d g0;
  q2g(g0, from.c);

  paths.clear();
  vector<Vector3d>::const_iterator it;
  for (it = vs.begin(); it != vs.end(); ++it) {    
    const Vector3d &v = *it;
    // reverse time if fwd=false
    GridPath<3> path;
    if (!Flow(path, g0, (fwd ? dt : -dt)*v))
      continue;
    
    paths.push_back(path);
  }
  return true;
}

