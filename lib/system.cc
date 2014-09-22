#include "system.h"
#include "state.h"
#include <assert.h>
#include <cmath>
#include <limits>
#include <cstring>

using namespace dsl;
using namespace std;

System::System(int n,
               int c)
{
  assert(n > 0);
  this->n = n;

  assert(c >= 0);
  this->c = c;

  this->xlb = new double[n];
  for (int i = 0; i < n; ++i)
    this->xlb[i] = std::numeric_limits<double>::min();

  this->xub = new double[n];
  for (int i = 0; i < n; ++i)
    this->xub[i] = std::numeric_limits<double>::max();

  this->ulb = new double[c];
  for (int i = 0; i < c; ++i)
    this->ulb[i] = std::numeric_limits<double>::min();

  this->uub = new double[c];
  for (int i = 0; i < c; ++i)
    this->uub[i] = std::numeric_limits<double>::max();
}


void System::SetBounds(const double *xlb,
                       const double *xub,
                       const double *ulb,
                       const double *uub)
{
  if (xlb)
    memcpy(this->xlb, xlb, n*sizeof(double));
  if (xub)
    memcpy(this->xub, xub, n*sizeof(double));
  if (ulb)
    memcpy(this->ulb, ulb, c*sizeof(double));
  if (uub)
    memcpy(this->uub, uub, c*sizeof(double));
}


System::~System()
{
  delete[] uub;
  delete[] ulb;
  delete[] xub;
  delete[] xlb;
}


State* System::Create(const double *x,
                      const double *u,
                      double t) const
{
  return new State(*this, x, u, t);
}


void System::Get(State &s, 
                 const State &sa, 
                 const State &sb, 
                 double a, 
                 bool time) const
{
  for (int j = 0; j < n; ++j)
    s.x[j] = (1 - a)*sa.x[j] + a*sb.x[j];  
  for (int j = 0; j < c; ++j)      
    s.u[j] = (1 - a)*sa.u[j] + a*sb.u[j];  
  if (time)
    s.t = (1 - a)*sa.t + a*sb.t;
}
