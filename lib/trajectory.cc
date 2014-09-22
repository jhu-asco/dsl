#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include <iomanip>
#include <cstring>
#include "trajectory.h"

using namespace dsl;
using namespace std;


Trajectory::Trajectory(const System &sys, int sn) : sys(sys)
{
  // the basic initialization always assumes that there is at least one state
  // normally only trim trajectories remain described by their initial state
  // since their dynamics is invariant and can be reconstructed directly

  this->sn = -1;
  this->states = 0;
  
  if (sn >= 0)
    this->Resize(sn);

  this->ext = false;  
  
  this->pn = 0;
  this->ps = 0;
  this->pmin = 0;
  this->pmax = 0;

  this->hmin = 1e-10;
}





namespace dsl {
  std::ostream& operator<<(std::ostream &os, const Trajectory &traj)
  {
    os << traj.sn << " " << traj.ext << endl;    
    for (int i = 0; i <= traj.sn; ++i)
      os << *traj.states[i] << endl;    
    return os;
  };
  
  std::istream& operator>>(std::istream &is, Trajectory &traj)
  {    
    if (traj.states) {
      for (int i = 0; i <= traj.sn; ++i)
        delete traj.states[i];
      free(traj.states);      
    }

    is >> traj.sn  >> traj.ext;
    traj.states = (State**)malloc((traj.sn+1)*sizeof(State*));    
    for (int i = 0; i <= traj.sn; ++i) {
      traj.states[i] = traj.sys.Create();
      is >> *traj.states[i];
    }
    return is;
  }
}


Trajectory::Trajectory(const Trajectory &traj) : sys(traj.sys)
{
  this->sn = traj.sn;
  this->states = (State**)malloc((traj.sn+1)*sizeof(State*));    
  for (int i = 0; i <= this->sn; ++i)
    this->states[i] = traj.states[i]->Clone();    
  
  this->ext = traj.ext;
  
  this->pn = traj.pn;
  this->ps = new double[traj.pn];
  this->pmin = new double[traj.pn];
  this->pmax = new double[traj.pn];
  if (traj.pn) {
    memcpy(this->ps, traj.ps, traj.pn*sizeof(double));
    memcpy(this->pmin, traj.pmin, traj.pn*sizeof(double));
    memcpy(this->pmax, traj.pmax, traj.pn*sizeof(double));
  }

  this->hmin = 1e-10;
}


Trajectory& Trajectory::operator=(const Trajectory &traj)
{
  if (&traj == this)
    return *this;
  
  assert(&traj.sys == &this->sys);
  
  if (this->sn != traj.sn) {
    // if not the same length then adjust states
    for (int i = 0; i <= this->sn; ++i)
      delete this->states[i];
    free(this->states);
    
    this->sn = traj.sn;

    this->states = (State**)malloc((this->sn+1)*sizeof(State*));
    //    this->states = new State*[this->sn + 1];
    for (int i = 0; i <= this->sn; ++i)
      this->states[i] = traj.states[i]->Clone();    
  } else {
    // otherwise just copy states
    for (int i = 0; i <= this->sn; ++i)
      *this->states[i] = *traj.states[i];        
  }

  this->ext = traj.ext;  
  
  if (this->pn != traj.pn) {
    delete[] this->ps;
    delete[] this->pmin;
    delete[] this->pmax;

    this->pn = traj.pn;
    
    this->ps = new double[traj.pn];
    this->pmin = new double[traj.pn];
    this->pmax = new double[traj.pn];
  }

  if (traj.pn) {
    memcpy(this->ps, traj.ps, traj.pn*sizeof(double));
    memcpy(this->pmin, traj.pmin, traj.pn*sizeof(double));
    memcpy(this->pmax, traj.pmax, traj.pn*sizeof(double));
  }
  return *this;
}


Trajectory::Trajectory(const System &sys, std::istream& istr) : sys(sys)
{
  istr >> sn  >> ext;
  states = (State**)malloc((sn+1)*sizeof(State*));    
  //  states = new State*[sn+1];
  
  for (int i=0; i<=sn; ++i) {
    states[i] = sys.Create();
    istr >> *states[i];
  }
  
  this->pn = 0;
  this->pmin = 0;
  this->pmax = 0;

  this->hmin = 1e-10;
}


Trajectory::~Trajectory()
{
  delete[] pmax;
  delete[] pmin;
  delete[] ps;

  for (int i=0; i<=sn; ++i)
    delete states[i];
  free(states);
}



Trajectory* Trajectory::Clone() const
{
  return new Trajectory(*this);
}


void Trajectory::Resize(int sn)
{
  if (sn == this->sn)
    return;

  for (int i = sn + 1; i <= this->sn; ++i)
    delete states[i];  
  
  states = (State**)realloc(states, (sn + 1)*sizeof(State*));

  for (int i = this->sn+1; i <= sn; ++i)
    states[i] = sys.Create();

  this->sn = sn;
}

void Trajectory::Reverse()
{
  State s(sys);
  for (int i = 0; i <= sn/2; ++i) {
    s = *states[i];
    *states[i] = *states[sn - i];
    *states[sn - i] = s;
  }
}




void Trajectory::Init(int sn,
                      const State *si,
                      const State *sf)
{
  assert(sn >= 0);
  Resize(sn);  
  if (si && sf) {
    for (int i = 0; i <= sn; ++i)      
      sys.Get(*states[i], *si, *sf, ((double)i)/sn);
  }
}


void Trajectory::Attach(const Trajectory &traj,
                        bool back,
                        bool time,
                        bool js)
{
  int sn = this->sn;        // old size
  
  if (sn < 0)
    js = 0;

  Resize(sn + traj.sn + 1 - js);
  
  //  memcpy(states + (sn + 1), traj.states + js, (traj.sn + 1 - js)*sizeof(State*));
  
  if (back) {
    double t = (sn < 0) ? 0 : states[sn]->t;
    int j = sn + 1;
    for (int i = js; i <= traj.sn; ++i) {
      *this->states[j] = *traj.states[i];     // copy states
      if (time)
        states[j]->t += t;                 // increment time
      ++j;
    }
    assert(j == this->sn + 1);
  } else {
    memmove(states + traj.sn + 1 - js, states, (sn + 1)*sizeof(State*));
    for (int i = 0; i <= traj.sn - js; ++i) {      
      *states[i] = *traj.states[i];     // copy states in the front
    }
    if (time) {
      double t = traj.states[traj.sn]->t;
      for (int i = traj.sn + js; i <= sn; ++i) {              
        states[i]->t += t;        // decrement time
      }
    }
  }
}


void Trajectory::Append(const State &s)
{
  Resize(sn + 1);
  *states[sn] = s;
}



void Trajectory::Clear()
{
  for (int i = 0; i <= sn; ++i)
    delete states[i];
  free(states);
  
  this->sn = -1;
  this->states = 0;
}


void Trajectory::Get(State &s) const
{
  if (!sn)
    return;
  
  // find the time slot
  if ( (s.t < states[0]->t - hmin) || (s.t > states[sn]->t + hmin)) {
    cerr << "[W] dgc::Trajectory::GetState:\t time t= " << s.t << " out of bounds [" << states[0]->t << "," << states[sn]->t << "]!" << endl;
    return;
  }


  // assume constant timestep  
  double h = (states[1]->t - states[0]->t);

  if (h < hmin) {
    cerr << "[W] dgc::Trajectory::GetState:\t time step is very small h=" << h << "!" << endl;    
  }


  if (fabs(s.t - states[0]->t) < hmin) {
    memcpy(s.x, states[0]->x, sys.n*sizeof(double));
    memcpy(s.u, states[0]->u, sys.c*sizeof(double));
    return;
  }
  
  if (fabs(s.t - states[sn]->t) < hmin) {
    memcpy(s.x, states[sn]->x, sys.n*sizeof(double));
    memcpy(s.u, states[sn]->u, sys.c*sizeof(double));
    return;
  }

  /*
  double i = (s.t - states[0]->t)/h;
  int i0 = MAX((int)floor(i), 0);
  int i1 = MIN((int)ceil(i), sn);
  cout << "t=" << s.t << " ti=" << states[0]->t << " tf=" << states[sn]->t << endl;
  cout << "i1=" << i1 << " sn=" << sn << " i0=" << i0 << endl;
  assert(i0 >= 0 && i0 <= sn);
  assert(i1 >= 0 && i1 <= sn);
  double d = i - i0;
  
  for (int j=0; j < sys.n; ++j)
    s.x[j] = (1-d)*states[i0]->x[j] + d*states[i1]->x[j];
  
  for (int j=0; j < sys.c; ++j)
    s.u[j] = (1-d)*states[i0]->u[j] + d*states[i1]->u[j];
  */

  double i = (s.t - states[0]->t)/h;
  int i0 = (int)i;

  //  cout << "t=" << s.t << " ti=" << states[0]->t << " tf=" << states[sn]->t << " h=" << h << endl;
  //  cout << " sn=" << sn << " i=" << i << " i0=" << i0 << endl;

  assert(i0 >= 0 && i0 < sn);

  //  cout << "i0=" << i0 << " i=" << i <<  " sn=" << sn << endl;
  
  sys.Get(s, *states[i0], *states[i0+1], i - i0, false);
}


int Trajectory::Get(double t) const
{
  // find the time slot
  if ( (t < states[0]->t - hmin) || (t > states[sn]->t + hmin)) {
    cerr << "[W] dgc::Trajectory::Get:\t time t= " << t << " out of bounds [" << states[0]->t << "," << states[sn]->t << "]!" << endl;
    return -1;
  }
  
  // assume constant timestep  
  double h = (states[1]->t - states[0]->t);
  
  if (h < hmin) {
    cerr << "[W] dgc::Trajectory::GetState:\t time step is very small h=" << h << "!" << endl;
  }
  
  if (fabs(t - states[0]->t) < hmin)
    return 0;
  
  if (fabs(t - states[sn]->t) < hmin)
    return sn-1;
  
  return (int)((t - states[0]->t)/h);
}



void Trajectory::Get(Trajectory &traj, double ti, double tf) const
{
  assert(traj.sn > 0);
  assert(tf > ti);
  double h = (tf - ti)/traj.sn;
  double t = ti;
  for (int i = 0; i <= traj.sn; ++i) {
    traj.states[i]->t = t;
    Get(*traj.states[i]);
    t += h;
  }
}



bool Trajectory::IsValidTime(double t) const
{
  if (t < states[0]->t || t > states[sn]->t)
    return false;
  return true;
}


void Trajectory::Resample(int sn)
{
  if (this->sn == sn)
    return;
  
  //  cout << "traj.sn=" << traj.sn << " dt=" << traj.states[traj.sn]->t - traj.states[0]->t << endl;

  State **states = (State**)malloc((sn+1)*sizeof(State*));
  //  State **states = new State*[sn + 1];
  double h = (this->states[this->sn]->t - this->states[0]->t)/sn;
  for (int i = 0; i <= sn; ++i) {
    states[i] = new State(this->sys);
    states[i]->t = this->states[0]->t + i*h;
    this->Get(*states[i]);
  }
  for (int i = 0; i <= this->sn; ++i)
    delete this->states[i];
  free(this->states);
  
  this->states = states;
  this->sn = sn;
}


void Trajectory::Refine(double s)
{
  int mi[this->sn];
  double mu[this->sn];
  for (int i = 0; i < this->sn; ++i) {
    mi[i] = i;
    mu[i] = s;
  }
  Modify(this->sn, mi, mu);
}


void Trajectory::Modify(int mn, const int *mi, const double *mu)
{
  assert(mi);
  assert(mu);

  State **states = (State**)malloc((sn + 1 + mn)*sizeof(State*));    
  //  State** states = new State*[this->sn + 1 + mn];
  
  int oj = 0;   // curent index in the original qs
  int nj = 0;   // current index in new qs
  for (int i = 0; i < mn; ++i) {
    assert(mi[i] < this->sn);    
    if (oj <= mi[i]) {
      memcpy(states + nj, this->states + oj, (mi[i] - oj + 1)*sizeof(State*));
      nj += (mi[i] - oj + 1);
      oj = mi[i] + 1;
    }
    assert(nj <= this->sn + mn);
    states[nj] = sys.Create();
    
    sys.Get(*states[nj], *this->states[oj-1], *this->states[oj], mu[i]);

    nj++;
    assert(nj <= this->sn + mn);
  }
  if (oj <= this->sn) {
    memcpy(states + nj, this->states + oj, (this->sn - oj + 1)*sizeof(State*));
  }
  
  this->sn += mn;
  free(this->states);
  this->states = states;
}



void Trajectory::Log(const char *logName, int di)
{
  std::ofstream file;
  file.open(logName);
  
  if (!file.is_open()) {
    std::cerr << "[E] dgc::Trajectory::Log: unable to open file" << logName << std::endl;
    return;
  }
  
  file.flags(ios::fixed);
  file.precision(8);  
  Log(file, di);
  file.close();
}


void Trajectory::Log(std::ofstream& os, int di)
{
  int fw = 20;
  for (int i = 0; i<=sn; i+=di) {
    os << setw(fw) << i;
    os << setw(fw) << states[i]->t;
    for (int j=0; j<sys.n; ++j)
      os << setw(fw) << states[i]->x[j];
   for (int j=0; j<sys.c; ++j)
      os << setw(fw) << states[i]->u[j];
    os << std::endl;
  }
}


void Trajectory::Log(std::ofstream& os, int fw, double a)
{
  if (isnan(a) || isinf(a) || fabs(a) > 1e8)
    a = 0; 
  os << setw(fw) << a;
}


void Trajectory::SetTime(double t0, double h)
{
  //  assert(h > 0);
  states[0]->t = t0;
  for (int i=1; i<=sn; ++i)
    states[i]->t = states[i-1]->t + h;
}

void Trajectory::SetParams(int pn, const double *ps, const double *pmin, const double *pmax)
{
  assert(pn > 0);
  this->pn = pn;
  this->ps = new double[pn];
  if (!ps) {
    for (int i = 0; i < pn; ++i)
      this->ps[i] = 0;
  } else {
    memcpy(this->ps, ps, pn*sizeof(double));
  }

  this->pmin = new double[pn];
  if (pmin) {
    memcpy(this->pmin, pmin, pn*sizeof(double));
  } else {
    for (int i = 0; i < pn; ++i)
      this->pmin[i] = 0;
  }
  this->pmax = new double[pn];
  if (pmax) {
    memcpy(this->pmax, pmax, pn*sizeof(double));
  } else {
    for (int i = 0; i < pn; ++i)
      this->pmax[i] = 1;
  }
}
