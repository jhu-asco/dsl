#include "carcost.h"
#include "utils.h"

namespace dsl {
using namespace std;

CarCost::CarCost(const CarGrid& grid, double ac, double eps)
:grid_(grid), use_twistnorm_metric_(false),wt_(Vector3d(0.1,1,2)), ac_(ac),eps_(eps){
  assert(eps > 0 && eps < 1);
  assert(ac >= 0);
}

CarCost::CarCost(const CarGrid& grid, const Vector3d& wt, double eps )
:grid_(grid), use_twistnorm_metric_(true), wt_(wt), ac_(0.1),eps_(eps){
  assert(eps > 0 && eps < 1);
  if(!(eps > 0 && eps < 1)){
    cout<<"wrong eps in carcost"<<endl;
    eps = 1e-6;
  }

  assert( (wt(0)>0) && (wt(1)>0) && (wt(2)>0) );
  if(!((wt(0)>0) && (wt(1)>0) && (wt(2)>0))){
    cout<<"wrong wt in carcost"<<endl;
    this->wt_ << 0.1, 1, 2;
  }
}

double CarCost::Real(const SE2Cell& a, const SE2Cell& b) const{
  //Find twist that takes you from a -> b
  const Matrix3d& ga = a.data;
  const Matrix3d& gb = b.data;
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg);

  //find cells along the way from a->b. If all the cell along the way don't exist then
  //it returns the max cost possible. Basically indicating they are not connected
  double d = twist.tail<2>().norm(); // total distance along curve
  int n_seg = ceil(d/(2*grid_.cs()(1))); // 2 * grid.cs[1] is to improve efficiency
  double d_seg = d/n_seg;
  SE2CellPtr wp; //waypoint cells
  for (int i_seg=1; i_seg<n_seg; i_seg++) {
    se2_exp(dg, (d_seg*i_seg / d) * twist);
    Matrix3d g = ga * dg;
    Vector3d axy; se2_g2q(axy, g);
    wp = grid_.Get(axy); //waypoint
    if (!wp)
      return numeric_limits<double>::quiet_NaN();
  }

  double wl = (twist.array()* wt_.array()).matrix().norm();//weighted length

  //Determine cost
  if(use_twistnorm_metric_){
    return wl;
  }else{
    return ac_*abs(twist(0)) + twist.tail<2>().norm();

  }
}

double CarCost::Heur(const SE2Cell& a, const SE2Cell& b) const {
  //Find twist that takes you from a -> b
  const Matrix3d& ga = a.data;
  const Matrix3d& gb = b.data;
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg);
  double wl = (twist.array()* wt_.array()).matrix().norm();//weighted length

  if(use_twistnorm_metric_){
    return (1 - eps_)*wl;
  }else{
    return (1-eps_)*(ac_*abs(twist(0)) + twist.tail<2>().norm());

  }
}
}
