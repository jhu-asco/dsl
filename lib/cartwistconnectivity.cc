#include <cartwistconnectivity.h>
#include <iostream>
#include "utils.h"
#include <algorithm>


namespace dsl {

using namespace Eigen;
using namespace std;

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Affine3d;
using std::vector;

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXbool;

CarTwistConnectivity::CarTwistConnectivity(const CarGrid& grid,const CarCost& cost, CarPrimitiveConfig& cfg)
:grid(grid),cost(cost){
  SetPrimitives(1, cfg);
}

CarTwistConnectivity::CarTwistConnectivity(const CarGrid& grid,const CarCost& cost,
                                 const vector<Eigen::Vector3d> &vs,
                                 double dt) : grid(grid), vs(vs), dt(dt),cost(cost)
{
}

CarTwistConnectivity::CarTwistConnectivity(const CarGrid& grid,const CarCost& cost,
                                 double dt,
                                 double vx,
                                 double kmax,
                                 int kseg,
                                 bool onlyfwd)
: grid(grid),cost(cost) {
  SetPrimitives(dt, vx, kmax, kseg, onlyfwd);
}

bool CarTwistConnectivity::SetPrimitives(double dt, CarPrimitiveConfig& cfg){

  cfg.nl = cfg.nl<1 ? 1: cfg.nl;  //! Make sure the number of different traj length is atleast 2
  cfg.na = cfg.na<3 ? 3: cfg.na;  //! Make sure the number of differentsteering angle is atleast 3

  double u = 1.0; //The speed can be anything. Using for clarity

  double tmin = cfg.lmin/u;
  double tmax = cfg.lmax/u;
  double wmax = u*cfg.tphioverlmax;
  double rmax = cfg.lmax/cfg.amax;
  double ymax = 2*(rmax - rmax*cos(cfg.amax));
  double xmax = 2*cfg.lmax;

  this->dt = dt;

  double del_t;
  if(cfg.nl==1)
    del_t = 1;
  else
    del_t = exp(log(tmax/tmin)/(cfg.nl-1));


  double del_w = wmax/(cfg.na-1);

  double sx = grid.cs(1);
  double sy = grid.cs(2);

  //create a grid which can accomodate any primitive starting from any orientation
  uint grid_nhc = ceil(max(ymax/sy, xmax/sx)); // number of half cells
  uint grid_nc  = 2*grid_nhc + 1;
  Matrix<bool, Dynamic, Dynamic> grid_seen(grid_nc, grid_nc);

  Affine3d igorg_to_mgcorg = Scaling(Vector3d(1/sx, 1/sy, 1)) *Translation3d(Vector3d(grid_nhc*sx, grid_nhc*sx, 0));

  //Allocate space for the
  vbs_per_angle.resize(grid.gs(0));
  for(size_t i=0;i<vbs_per_angle.size();i++)
    vbs_per_angle[i].reserve(cfg.na*cfg.nl);

  if(cfg.tocenter){
    for(int idx_a=0; idx_a < grid.gs(0);idx_a++){
      double th = grid.CellCenterIth(idx_a,0);
      for (size_t i = 0, size = grid_seen.size(); i < size; i++)
        *(grid_seen.data()+i) = false;
      Affine3d igorg_to_car = igorg_to_mgcorg * AngleAxisd(th,Vector3d::UnitZ());
      for(int idx_t=0 ; idx_t< cfg.nl;idx_t++){
        double t = tmin*pow(del_t, idx_t);
        for(int idx_w=0 ; idx_w < cfg.na; idx_w++){
          double w = idx_w*del_w;
          double tpert;
          if(cfg.pert)
            tpert = t*(1 + 0.2*(cos(40*w)-1)); // t perturbed
          else
            tpert = t;
          double apert = w*tpert;              // angle perturbed
          double lpert = u*tpert;
          if(apert > cfg.amax || lpert > 1.1*cfg.lmax || lpert < 0.7*cfg.lmin)
            continue;

          Matrix3d gend; se2_exp(gend,Vector3d(w*tpert, u*tpert, 0));
          Vector3d xyzend(gend(0,2),gend(1,2),0);
          Vector3d idxd = igorg_to_car * xyzend;
          Vector3i idxi(round(idxd(0)),round(idxd(1)),round(idxd(2)));

          if((uint)idxi(1)>grid_nc-1 || (uint)idxi(0)>grid_nc-1) //Shouldn't be necessary
            continue;

          if(grid_seen(idxi(1), idxi(0))==true)
            continue;
          else
            grid_seen(idxi(1), idxi(0))=true;

          Vector3d xyzend_snapped = igorg_to_car.inverse()* (idxi.cast<double>()) ;
          //cout<<"xyzend_snapped:"<<xyzend_snapped.transpose()<<endl;
          Vector2d WVx = getWVx(xyzend_snapped(0), xyzend_snapped(1));
          Vector3d vfp(  WVx(0),  WVx(1),0);
          Vector3d vfn( -WVx(0),  WVx(1),0);
          Vector3d vbp(  WVx(0), -WVx(1),0);
          Vector3d vbn( -WVx(0), -WVx(1),0);

          if(abs(WVx(0))<1e-10){
            vbs_per_angle[idx_a].push_back(vfp);
            if(!cfg.fwdonly)
              vbs_per_angle[idx_a].push_back(vbp);
          }else{
            vbs_per_angle[idx_a].push_back(vfp);
            vbs_per_angle[idx_a].push_back(vfn);
            if(!cfg.fwdonly){
              vbs_per_angle[idx_a].push_back(vbp);
              vbs_per_angle[idx_a].push_back(vbn);
            }
          }
        }
      }
    }
  }else{

    for(int idx_t=0 ; idx_t< cfg.nl;idx_t++){
      double t = tmin*pow(del_t, idx_t);
      for(int idx_w=0 ; idx_w < cfg.na; idx_w++){
        double w = idx_w*del_w;
        double tpert;
        if(cfg.pert)
          tpert = t*(1 + 0.2*(cos(40*w)-1)); // t perturbed
        else
          tpert = t;
        double apert = w*tpert;              // angle perturbed
        double lpert = u*tpert;
        if(apert > cfg.amax || lpert > cfg.lmax || lpert < 0.7*cfg.lmin)
          continue;

        Vector3d vfp(w*tpert, u*tpert,0);
        Vector3d vfn(-w*tpert, u*tpert,0);
        Vector3d vbp(w*tpert, -u*tpert,0);
        Vector3d vbn(-w*tpert, -u*tpert,0);
        //Put the same set of primitive for all angles
        for(int idx_a=0; idx_a < grid.gs(0);idx_a++){
          if(abs(w)<1e-10){
            vbs_per_angle[idx_a].push_back(vfp);
            if(!cfg.fwdonly)
              vbs_per_angle[idx_a].push_back(vbp);
          }else{
            vbs_per_angle[idx_a].push_back(vfp);
            vbs_per_angle[idx_a].push_back(vfn);
            if(!cfg.fwdonly){
              vbs_per_angle[idx_a].push_back(vbp);
              vbs_per_angle[idx_a].push_back(vbn);
            }
          }
        }
      }
    }

  }

  return true;
}

bool CarTwistConnectivity::SetPrimitives(double dt, double vx, double kmax, int kseg, bool onlyfwd) {
  if (dt <= 0)
    return false;

  kseg = kseg < 1 ? 1 : kseg;

  this->dt = dt;

  vs.clear();

  for (int i = 0; i <= kseg; i++) {
    double k = i*kmax/kseg;
    double w = vx*k;
    vs.push_back(Vector3d(w, vx, 0));
    vs.push_back(Vector3d(-w, vx, 0));
    if (!onlyfwd) {
      vs.push_back(Vector3d(w, -vx, 0));
      vs.push_back(Vector3d(-w, -vx, 0));
    }
  }

  return true;
}

bool CarTwistConnectivity::
operator()(const SE2Cell& from,
           std::vector< std::tuple<SE2CellPtr, SE2Twist, double> >& paths,
           bool fwd) const {
  const vector<Vector3d>* p_vstemp;
  int idx_a = grid.Index(from.c,0);
  if(!vbs_per_angle.empty())
    p_vstemp = &(vbs_per_angle[idx_a]);
  else
    p_vstemp = &vs;

  paths.clear();

  for (auto& s : *p_vstemp) {
    Vector3d v = s*(fwd ? dt: -dt);
    Matrix3d dg; se2_exp(dg, v);
    Matrix3d g_to; g_to = from.data*dg;
    Vector3d axy_to; se2_g2q(axy_to, g_to);
    SE2CellPtr to; to = grid.Get(axy_to);
    if (!to) // "to" cell is not a part of the grid(because it has an obstacle)
      continue;
    g_to = to->data;

    double primcost = cost.Real(from,*to);
    if(std::isnan(primcost)) //path is not clear
      continue;

    //no slip primitive that takes you from a cell to successors cell(only position not orientation)
    Matrix3d g_inv;
    if(fwd){
      se2_inv(g_inv,from.data);
      dg = g_inv*g_to;
    }else{
      se2_inv(g_inv,g_to);
      dg = g_inv*from.data;
    }
    Vector3d v_noslip;
    v_noslip << getWVx(dg(0,2),dg(1,2)),0; //no slip velocity to successors cell(only position not orientation)

    std::tuple< SE2CellPtr, SE2Twist, double> pathTuple(to,v_noslip,primcost);
    paths.push_back(pathTuple);
  }
  return true;
}


bool CarTwistConnectivity::GetPrims(const Vector3d pos, vector<vector<Vector2d>>& prims ){
  //Display the primitive at start
  SE2CellPtr cell_start = grid.Get(pos);
  vector<std::tuple< SE2CellPtr, SE2Twist, double> > paths;
  if(cell_start){
    (*this)(*cell_start,paths,true);
    Matrix3d g0 = cell_start->data;
    prims.reserve(paths.size());
    for(auto& path:paths){
      SE2CellPtr cell_to = std::get<0>(path);
      if(!cell_to)
        continue;

      SE2Twist& v=  std::get<1>(path);
      double d = fabs(v[1]); // total distance along curve
      int n_seg = ceil(d/ (2 * grid.cs[1])); // 2 * grid.cs[1] is to improve efficiency
      double s = d/n_seg;
      vector<Vector2d> prim(0); prim.reserve(n_seg+1);
      prim.push_back(Vector2d(g0(0,2),g0(1,2)));
      for (int i_seg=1; i_seg<=n_seg; i_seg++) {
        Vector3d axy;
        Matrix3d g, dg;
        se2_exp(dg, (s*i_seg / d) * v);
        g = g0 * dg;
        prim.push_back(Vector2d(g(0,2),g(1,2)));
      }
      prims.push_back(prim);
    }
    return true;
  }else{
    return false;
  }
}

}
