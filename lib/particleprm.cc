#include "particleprm.h"
#include <limits>
#include <iostream>
#include <math.h>
#include <assert.h>

using namespace dsl;
using namespace std;


#ifndef MAX
#define MAX(a,b)((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b)((a)<(b)?(a):(b))
#endif

#define DIST(a,b) (sqrt(((a)[0]-(b)[0])*((a)[0]-(b)[0])+((a)[1]-(b)[1])*((a)[1]-(b)[1])))
#define NORM(a) (sqrt((a)[0]*(a)[0]+(a)[1]*(a)[1]))
#define NORM3(a) (sqrt((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2]))


ParticlePrm::ParticlePrm(const State &slb,
                         const State &sub) : 
  ds(1), eps(1e-6), tf(0), h(.1)
{
}

ParticlePrm::~ParticlePrm()
{
}


bool ParticlePrm::sint_opt(double ts[2], double as[2], 
                           double xi, double vi, 
                           double xf, double vf,
                           double *am, int n, double tf)
{

  // solution for t2
  //-(2*vf + 2^(1/2)*(vf^2 + vi^2 + 2*a*xf - 2*a*xi)^(1/2))/(2*a)
  //-(2*vf - 2^(1/2)*(vf^2 + vi^2 + 2*a*xf - 2*a*xi)^(1/2))/(2*a)
  // t1 = t2+(vf-vi)/a

  //cout << xi << " " << xf << " " << vi << " " << vf  << " " << n << endl;
  
  double tmin = numeric_limits<double>::max();

  double a, cp, t1, t2;

  for (int i = 0; i < n; ++i) {
    a = am[i];
  
    cp = 2*(vf*vf + vi*vi + 2*a*xf - 2*a*xi);
    //cout << cp <<  " " << a << endl;
    if (cp >= 0) {
      t2 =  -(2*vf + sqrt(cp))/(2*a);
      //cout << "t2: " << t2 << endl;
      if (t2 > 0) {
        t1 = t2 + (vf-vi)/a;
        //cout << "t1: " << t1 << endl;
        if (t1 > 0) {
          double d = fabs(t1 + t2 - tf);
          //cout << "d: " << d << endl;
          if (d < tmin) {
            tmin = d;
            ts[0] = t1; ts[1] = t2;
            as[0] = a; as[1] = -a;
          }
        }
      }
      
      t2 =  -(2*vf - sqrt(cp))/(2*a);
      //cout << "t2: " << t2 << endl;
      if (t2 > 0) {
        t1 = t2 + (vf-vi)/a;
        //cout << "t1: " << t1 << endl;
        if (t1 > 0) {
          double d = fabs(t1 + t2 - tf);
          //cout << "d: " << d << endl;
          if (d < tmin) {
            tmin = d;
            ts[0] = t1; ts[1] = t2;
            as[0] = a; as[1] = -a;
          }
        }    
      }
    }
  }

  if (tmin == numeric_limits<double>::max())
    return false;

  return true;
}


bool ParticlePrm::sint_a(double a[2], 
                         double xi, double vi, 
                         double xf, double vf, double tf)
{
  double c = 2*tf*tf*vf*vf + 2*tf*tf*vi*vi - 4*tf*vf*xf + 4*tf*vf*xi - 4*tf*vi*xf + 4*tf*vi*xi + 4*xf*xf - 8*xf*xi + 4*xi*xi;
  
  if (c >= 0) {
    a[0] = (2*xf - 2*xi - tf*(vf + vi) - sqrt(c))/(tf*tf);
    a[1] = (2*xf - 2*xi - tf*(vf + vi) + sqrt(c))/(tf*tf);
    return true;
  } else {
    return false;
  }
}

double ParticlePrm::sint_profile(double tss[][2], double ass[][2], 
                                 const State &sa, const State &sb)
{
  int d = 3;

  double imax = 0;
  double tmax = 0;
  //  double tss[d][2];
  //  double ass[d][2];
  
  double ams[2] = {sa.sys.uub[0], -sa.sys.uub[0]};

  double *xi = sa.x;
  double *vi = sa.x+d;
  double *xf = sb.x;
  double *vf = sb.x+d;

  for (int i = 0; i < d; ++i) {
    if (fabs(xi[i]-xf[i]) + fabs(vi[i]-vf[i]) < 1e-6) {
      tss[i][0] = tss[i][1] = 0;
      ass[i][0] = ass[i][1] = 0;
      continue;
    }
    
    //assert(sint_opt(tss[i], ass[i], xi[i], vi[i], xf[i], vf[i], ams, 2, 0));
    if(!sint_opt(tss[i], ass[i], xi[i], vi[i], xf[i], vf[i], ams, 2, 0))
      return -1;
    double tt = tss[i][0] + tss[i][1];
    if (tt > tmax) {
      tmax = tt;
      imax = i;
    }
  }
  
  double am[2];
  for (int i = 0; i < d; ++i) {    
    if (i != imax) {
      if (fabs(xi[i]-xf[i]) + fabs(vi[i]-vf[i]) < 1e-6) {
        tss[i][0] = tss[i][1] = 0;
        ass[i][0] = ass[i][1] = 0;
        continue;
      }
      //assert(sint_a(am, xi[i], vi[i], xf[i], vf[i], tmax));
      if(!sint_a(am, xi[i], vi[i], xf[i], vf[i], tmax))
        return -1;
      
      bool ok = sint_opt(tss[i], ass[i], xi[i], vi[i], xf[i], vf[i], am, 2, tmax);
      if (!ok) {
        //cout << "[W] ParticlePrm::Compute: could not solve sint_opt b/n:";
        //cout << "sa=" << sa << endl;
        //cout << "sb=" << sb << endl;
        return -1;
      }
    }
  }
  return tmax;
}


bool ParticlePrm::Compute(Trajectory &traj, const State &sa, const State &sb)
{
  //  assert(NORM(sb.x+2) > 1e-10);
  
  //  if (sb.t < 0)
    //    cout << "sa.t=" <<  sa.t << " sb.t=" << sb.t << endl;

  if (this->tf > 0 && sb.t > this->tf)
    return false;

  int d = 3;

  const double *qa = sa.x;
  const double *qb = sb.x;
  const double *va = sa.x+d;
  const double *vb = sb.x+d;
  

  // shortest distance
  double l = DIST(qb, qa);
  if (l < eps)
    return false;
  
  // max speed
  double vmax = ( d == 2 ? NORM(traj.sys.xub + d) : NORM3(traj.sys.xub + d));

  //  cout << sa << "->" << sb << endl;  

  // check if final time is feasible
  if (sb.t > eps && sb.t - sa.t < l/vmax)
    return false;

  //  cout << "sb.t=" << sb.t << endl;

  bool kin = (traj.sys.uub[0] > 1/eps);

  if (kin) {

    bool ok = (sb.t < 0);

    //    double tb = sb.t;
    double tb = 0;
    //    cout << "TB=" << tb << endl;

    if (tb < eps)
      tb = sa.t + l/vmax;
    
    if (this->tf > 0 && tb > this->tf) {      
      return false;
    }

    double T = tb - sa.t;

    if (T < 2*this->h)
      return false;

    
    int kn = MAX(T/this->h, 1);
    traj.Resize(kn);
    
    double h = T/kn;
    
    //    if (ok)
    //      cout << "h=" << h << " sa=" << sa << " tb=" << sb.t << endl;

    assert(h > 0);
    
    *traj.states[0] = sa;
    for (int k = 1; k <= kn; ++k) {
      State &sp = *traj.states[k-1];
      State &s = *traj.states[k];
      s.t = sp.t + h;
      for (int i = 0; i < d; ++i) {
        s.x[i] = qa[i] + (s.t-sa.t)/T*(qb[i] - qa[i]);        
        sp.x[d + i] = (s.x[i] - sp.x[i])/h;
      }
    }

    //    ((State&)sb).t = traj.states[kn]->t;
    //    ((State&)sb).t = tb;
    
    return true;
  } 

  // time optimal double integrator
  
  double tss[d][2];
  double ass[d][2];

  double tmax = sint_profile(tss, ass, sa, sb);

  //  if (tmax < 0)
  //    return false;
  //  assert(tmax > 0);
  
  if (tmax < 2*this->h)
    return false;
  
  int kn = MAX(tmax/this->h, 1);
  traj.Resize(kn);
    
  double h = tmax/kn;
  
  *traj.states[0] = sa;
  for (int k = 1; k <= kn; ++k) {
    State &s = *traj.states[k];
    double *q = s.x;
    double *v = s.x + d;
    s.t = k*h;

    for (int i = 0; i < d; ++i) {            

      const double &t1 = tss[i][0];
      const double &a1 = ass[i][0];
      double ts = s.t - t1;

      // new acc begins
      if (ts > 0) { // next segment
        double &a2 = ass[i][1];
        v[i] = va[i] + t1*a1 + ts*a2;
        q[i]  = qa[i] + va[i]*t1 + t1*t1*a1/2 + (va[i] + t1*a1)*ts + ts*ts*a2/2;
        s.u[i] = a2;              
      } else {
        v[i] = va[i] + a1*s.t;
        q[i]  = qa[i] + va[i]*s.t + s.t*s.t*a1/2;
        s.u[i] = a1;
      }
    }
  }

  //  cout << "ERR=" << DIST(traj.states[traj.sn]->x, sb.x) << endl;
  // kill numerical inaccuracies

  //  memcpy(traj.states[traj.sn]->x, sb.x, traj.sys.n*sizeof(double));
  
  return true;  

  assert(false);

  // spline

  if (sb.t < eps)
    ((State&)sb).t = sa.t + l/vmax;

  if (this->tf > 0 && sb.t > this->tf)
    return false;

  double c2[d];
  double c3[d];

  double T = sb.t - sa.t;
  
  if (T < 2*this->h)
    return false;

  kn = MAX(T/this->h, 1);
  traj.Resize(kn);

  h = T/kn;

  for (int i = 0; i < d; ++i) {
    c2[i] = 3*(qb[i] - qa[i])/(T*T) - (2*va[i] + vb[i])/T;
    c3[i] = -2*(qb[i] - qa[i])/(T*T*T) + (va[i] + vb[i])/(T*T);
  }

  double t = 0;
  for (int k = 0; k <= kn; ++k) {
    for (int i = 0; i < d; ++i) {
      traj.states[k]->x[i] = qa[i] + t*va[i] + t*t*c2[i] + t*t*t*c3[i];
      traj.states[k]->x[d + i] = va[i] + 2*t*c2[i] + 3*t*t*c3[i];
    }
    traj.states[k]->t = sa.t + t;
    t += h;
  }
  ((State&)sb).t = traj.states[kn]->t;

  //  cout << "sb=" << sb << endl;
  //  cout << "s_=" << *traj.states[traj.sn] << endl;

  /*
  double T = 1;;

  for (int i = 0; i < d; ++i) {
    c2[i] = 3*(qb[i] - qa[i])/(T*T) - (2*va[i] + vb[i])/T;
    c3[i] = -2*(qb[i] - qa[i])/(T*T*T) + (va[i] + vb[i])/(T*T);
  }

  //  int kn = MAX(1,(int)(DIST(qb,qa)/.2));

  int kn = MAX(l/.1, 1);
  traj.Init(kn);

  double h = T/kn;
  
  for (int k = 0; k <= kn; ++k) {
    double t = k*h;      
    for (int i = 0; i < d; ++i) {
      traj.states[k]->x[i] = qa[i] + t*va[i] + t*t*c2[i] + t*t*t*c3[i];
      traj.states[k]->x[d + i] = va[i] + 2*t*c2[i] + 3*t*t*c3[i];
    }
    if (!k) {
      traj.states[k]->t = ta;
    } else {
      const double *xa = traj.states[k - 1]->x;
      const double *xb = traj.states[k]->x;
      double ds = 0;
      double v = 0;
      for (int i = 0; i < d; ++i) {
        double dxi = xb[i] - xa[i];
        ds += (dxi*dxi);
        double vi = (xa[d + i] + xb[d + i])/2;
        v += vi*vi;
      }
      ds = sqrt(ds);
      v = sqrt(v);
      traj.states[k]->t = traj.states[k - 1]->t + ds/v;
    }      
  }

  ((State&)sb).t = traj.states[kn]->t;
  */

  return true;
}
