#include "state.h"
#include "system.h"
#include <string.h>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <stdlib.h>

using namespace dsl;
using namespace std;


State::State(const System& sys,
             const double *x,
             const double *u,
             double t) : 
  sys(sys)
{
  this->x = (double*)malloc(sys.n*sizeof(double));
  this->u = (double*)malloc(sys.c*sizeof(double));
  
  if (x)
    memcpy(this->x, x, sys.n*sizeof(double));
  else
    memset(this->x, 0, sys.n*sizeof(double));
  
  if (u)
    memcpy(this->u, u, sys.c*sizeof(double));
  else
    memset(this->u, 0, sys.c*sizeof(double));
  
  this->t = t;

}


State::State(const State &s) : sys(s.sys)
{

  this->x = (double*)malloc(sys.n*sizeof(double));
  this->u = (double*)malloc(sys.c*sizeof(double));
  
  (*this) = s;
}

State& State::operator=(const State &s)
{
  if (&s == this)
    return *this;

  assert(&this->sys == &s.sys);
      
  memcpy(this->x, s.x, sys.n*sizeof(double));  
  memcpy(this->u, s.u, sys.c*sizeof(double));
  this->t = s.t;

  return *this;
} 


State::State(const System& sys, std::istream& istr) :
  sys(sys)
{  
  this->x = (double*)malloc(sys.n*sizeof(double));
  this->u = (double*)malloc(sys.c*sizeof(double));

  istr >> t;
  
  for (int i = 0; i < sys.n; ++i)
    istr >> x[i];
  
  for (int i = 0; i < sys.c; ++i)
    istr >> u[i];
  
}

State* State::Clone() const
{
  return new State(*this);
}


State::~State()
{
  free(u);
  free(x);
}


namespace dsl {
  std::ostream& operator<<(std::ostream &os, const State &s)
  {
    os.precision(20);
    os << s.t;
    for (int i = 0; i < s.sys.n; ++i)
      os << std::setw(40) << s.x[i];
      //      os << std::scientific << std::setw(40) << s.x[i];
    for (int i = 0; i < s.sys.c; ++i)
      os << std::setw(40) << s.u[i];
      //      os << std::scientific << std::setw(40) << s.u[i];
    return os;
  };
  
  std::istream& operator>>(std::istream &is, State &s)
  {    
    is >> s.t;    
    for (int i = 0; i < s.sys.n; ++i)
      is >> s.x[i];    
    for (int i = 0; i < s.sys.c; ++i)
      is >> s.u[i];
    return is;
  }
}
