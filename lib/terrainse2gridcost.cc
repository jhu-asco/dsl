#include <terrainse2gridcost.h>
#include "utils.h"
#include <exception>

namespace dsl {
using namespace std;
using namespace Eigen;
TerrainSE2GridCost::TerrainSE2GridCost(const TerrainSE2Grid& grid, const SE2GridCostConfig& config)
  :grid_(grid),config_(config), trav_min_(numeric_limits<double>::max()){

  double& eps = config_.eps;
  Vector3d& wt = config_.twist_weight;

  assert(eps > 0 && eps < 1);
  if(!(eps > 0 && eps < 1)){
    throw invalid_argument( "Wrong eps in TerrainSE2GridCost. 0 < eps < 1 desired." );
    eps = 1e-6;
  }

  assert( (wt(0)>0) && (wt(1)>0) && (wt(2)>0) );
  if(!((wt(0)>0) && (wt(1)>0) && (wt(2)>0))){
    throw invalid_argument( "Wrong wt in TerrainSE2GridCost. wt(i) > 0 for all i." );
    wt << 0.1, 1, 2;
  }

  for(int id=0; id < grid.nc(); id++)
    if(grid.cells()[id])
      trav_min_ = trav_min_ < grid.cells()[id]->data.traversibility ?
                  trav_min_ : grid.cells()[id]->data.traversibility;

}

double TerrainSE2GridCost::Real(const SE2TerrainCell& a, const SE2TerrainCell& b) const {
  const Vector3d& wt = config_.twist_weight;
  const bool& sp = config_.subpixel;

  //Find twist that takes you from a -> b
  Matrix3d ga; se2_q2g(ga,a.c);
  Matrix3d gb; se2_q2g(gb,b.c);
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg);

  //find cells along the way from a->b. If all the cell along the way don't exist then
  //it returns nan, indicating they are not connected
  double d = twist.tail<2>().norm(); // total distance along curve
  int n_seg = ceil(d/ grid_.cs()[1]); // want segments to lie in each cell along the way
  double d_seg = d/n_seg;
  SE2TerrainCell::Ptr wp; //waypoint cells
  double trav_sum = a.data.traversibility + b.data.traversibility;
  double delh_pve_sum = 0; //positive change in height. Only positive change in height incurs cost.
  double pot_prev = a.data.height;
  for (int i_seg=1; i_seg<n_seg; i_seg++) {
    se2_exp(dg, (d_seg*i_seg / d) * twist);
    Matrix3d g = ga * dg;
    Vector3d axy; se2_g2q(axy, g);
    wp = grid_.Get(axy);
    if (!wp) //one of the waypoints is outside the grid.
      return numeric_limits<double>::quiet_NaN();

    trav_sum += wp->data.traversibility;
    if( (wp->data.height - pot_prev) > 0 )
      delh_pve_sum += wp->data.height - pot_prev;

    pot_prev = wp->data.height;
  }
  if( (b.data.height - pot_prev) > 0 )
    delh_pve_sum += b.data.height - pot_prev;

  double trav_avg = trav_sum/(n_seg+1);

  double wl = (twist.array()* wt.array()).matrix().norm();//weighted length

  return delh_pve_sum + trav_avg*wl;
}

double TerrainSE2GridCost::Heur(const SE2TerrainCell& a, const SE2TerrainCell& b) const {
  const double& eps = config_.eps;
  const Vector3d& wt = config_.twist_weight;

  Matrix3d ga; se2_q2g(ga,a.c);
  Matrix3d gb; se2_q2g(gb,b.c);
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg); //se2 twist that takes you from a -> b
  double wl = (twist.array()* wt.array()).matrix().norm();//weighted length
  double delh_pve_sum = (b.data.height - a.data.height) > 0 ?
                        (b.data.height - a.data.height) : 0; // only positive change in height incurs cost

  return (1-eps)*(delh_pve_sum + trav_min_ * wl);
}

}
