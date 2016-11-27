#include <terrainse2gridcost.h>
#include "utils.h"

namespace dsl {
TerrainSE2GridCost::TerrainSE2GridCost(const TerrainSE2Grid& grid, const Vector3d& w, double eps)
  :grid(grid),wt(wt),eps(eps){
  assert(eps > 0 && eps < 1);
  if(!(eps > 0 && eps < 1)){
    cout<<"wrong eps in carcost"<<endl;
    eps = 1e-6;
  }

  assert( (wt(0)>0) && (wt(1)>0) && (wt(2)>0) );
  if(!((wt(0)>0) && (wt(1)>0) && (wt(2)>0))){
    cout<<"wrong wt in carcost"<<endl;
    this->wt << 0.1, 1, 2;
  }
}

double TerrainSE2GridCost::Real(const TerrainCell& a, const TerrainCell& b) const {
  //Find twist that takes you from a -> b
  Matrix3d ga; se2_q2g(ga,a.c);
  Matrix3d gb; se2_q2g(gb,b.c);
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg);

  //find cells along the way from a->b. If all the cell along the way don't exist then
  //it returns nan, indicating they are not connected
  double d = twist.tail<2>().norm(); // total distance along curve
  int n_seg = ceil(d/ grid.cs[1]); // 2 * grid.cs[1] is to improve efficiency
  double d_seg = d/n_seg;
  TerrainCellPtr wp; //waypoint cells
  double trav_sum = a.data.traversibility + b.data.traversibility;
  for (int i_seg=1; i_seg<n_seg; i_seg++) {
    se2_exp(dg, (d_seg*i_seg / d) * twist);
    Matrix3d g = ga * dg;
    Vector3d axy; se2_g2q(axy, g);
    wp = grid.Get(axy);
    if (!wp)
      return numeric_limits<double>::max();
    trav_sum += wp->data.traversibility;
  }
  double trav_avg = trav_sum/(n_seg+1);

  double wl = (twist.array()* wt.array()).matrix().norm();//weighted length
  double pot = b.data.height - a.data.height; //change in potential energy

  return pot + trav_avg*wl;

}

double TerrainSE2GridCost::Heur(const TerrainCell& a, const TerrainCell& b) const {
  Matrix3d ga; se2_q2g(ga,a.c);
  Matrix3d gb; se2_q2g(gb,b.c);
  Matrix3d gai; se2_inv(gai,ga);
  Matrix3d dg = gai*gb;
  Vector3d twist; se2_log(twist,dg); //se2 twist that takes you from a -> b
  double wl = (twist.array()* wt.array()).matrix().norm();//weighted length
  double pot = b.data.height - a.data.height; //change in potential energy

  return pot+wl; //True cost is pot + traversibility*wl
}

}
