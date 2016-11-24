#include "carconnectivity.h"

namespace dsl{

template<>
bool CarTwistConnectivity::operator()(const SE2Cell& from,
                                           std::vector< std::tuple<SE2CellPtr, SE2Twist, double> >& paths,
                                           bool fwd) const{
  const vector<Vector3d>* vbstemp;
  int idx_a = grid.Index(from.c,0);
  if(!vbs_per_angle.empty())
    vbstemp = &(vbs_per_angle[idx_a]);
  else
    vbstemp = &vbs;

  paths.clear();

  for (auto& s : *vbstemp) {
    Vector3d vdt_suggested = s*(fwd ? dt: -dt); //suggested twist(because it will be snapped to cell center)
    Matrix3d dg; se2_exp(dg, vdt_suggested);//relative transform
    Matrix3d g_to; g_to = from.data*dg;//end pose resulting from the twist
    Vector3d axy; se2_g2q(axy, g_to);
    SE2CellPtr to; to = grid.Get(axy);
    if (!to) // "to" cell is not a part of the grid(because it has an obstacle)
      continue;

    //Get cost and check if path is clear
    double primcost = cost.Real(from,*to);
    if(std::isnan(primcost)) //path is not clear
      continue;

    //no slip primitive that takes you from a cell to successors cell(only position not orientation)
    Matrix3d g_inv;
    if(fwd){
      se2_inv(g_inv,from.data);
      dg = g_inv*to->data;
    }else{
      se2_inv(g_inv,to->data);
      dg = g_inv*from.data;
    }
    Vector3d v_noslip;
    v_noslip << getWVx(dg(0,2),dg(1,2)),0; //no slip velocity to successors cell(only position not orientation)

    std::tuple< SE2CellPtr, SE2Twist, double> pathTuple(to,v_noslip,primcost);
    paths.push_back(pathTuple);
  }
  return true;
}

template<>
bool CarPathConnectivity::
    operator()(const SE2Cell& from,
               std::vector< std::tuple<SE2CellPtr, SE2Path, double> >& paths,
               bool fwd) const {
  paths.clear();
  for (auto& vb : vbs) {
    Vector3d vdt_suggested = vb*(fwd ? dt: -dt); //suggested twist(because it will be snapped to cell center)
    Matrix3d dg; se2_exp(dg, vdt_suggested); //relative transform
    Matrix3d g; g = from.data*dg; //end pose resulting from the twist
    Vector3d axy; se2_g2q(axy, g);
    SE2CellPtr to; to = grid.Get(axy); //snap the end pose to grid
    if (!to) // "to" cell is not a part of the grid(because it has an obstacle)
      continue;

    //Get cost and check if path is clear
    double primcost = cost.Real(from,*to);
    if(std::isnan(primcost)) //path is not clear
      continue;

    //sample states along the path
    Matrix3d gi; se2_inv(gi,from.data); dg = gi*to->data; //relative of from.data to to->data
    Vector3d vdt_slip; se2_log(vdt_slip,dg);//twist that take you exactly to successor
    Vector3d vdt_noslip;vdt_noslip << getWVx(dg(0,2),dg(1,2)),0;//twist that takes you to successor(exactly only in position, not orientation))
    Vector3d vdt = vdt_slip;//chose between slip and no slip version
    double d = vdt.tail<2>().norm(); // total distance along curve
    int n_seg = ceil(d/(2*grid.cs[1])); // 2 * grid.cs[1] is to improve efficiency
    double d_seg = d/n_seg;
    SE2Path path; path.resize(n_seg); // path doesn't include start
    if(fwd){
      for (int i_seg=1; i_seg<=n_seg; i_seg++) {
        se2_exp(dg, (d_seg*i_seg / d) * vdt);
        path[i_seg-1] = from.data * dg;
      }
    }else{
      for (int i_seg=1; i_seg<=n_seg; i_seg++) {
        se2_exp(dg, -(d_seg*i_seg / d) * vdt);
        path[i_seg-1] = to->data * dg;
      }
    }

    //add primitive(along with cost and to-cell) to the paths
    std::tuple< SE2CellPtr, SE2Path, double> pathTuple(to,path,primcost);
    paths.push_back(pathTuple);
  }
  return true;
}
}
