#include "carconnectivity.h"
#include <iostream>

using namespace dsl;
using namespace std;


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
  
  v = 1;  // fixed forward velocity
  w = tphi*v; 

  vs.push_back(Vector3d(v, 0, 0));
  vs.push_back(Vector3d(v, 0, w));
  vs.push_back(Vector3d(v, 0, -w));
  vs.push_back(Vector3d(-v, 0, 0));
  vs.push_back(Vector3d(-v, 0, w));
  vs.push_back(Vector3d(-v, 0, -w));

  dt = 1;
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

