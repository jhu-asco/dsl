#include <carconntwistprim.h>
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

CarConnTwistPrim::CarConnTwistPrim(const CarGrid& grid, CarPrimitiveCfg& cfg)
:grid(grid){
  SetPrimitives(1, cfg);
}

CarConnTwistPrim::CarConnTwistPrim(const CarGrid& grid,
                                 const vector<Eigen::Vector3d> &vs,
                                 double dt) : grid(grid), vs(vs), dt(dt)
{
}

CarConnTwistPrim::CarConnTwistPrim(const CarGrid& grid,
                                 double dt,
                                 double vx,
                                 double kmax,
                                 int kseg,
                                 bool onlyfwd)
: grid(grid) {
  SetPrimitives(dt, vx, kmax, kseg, onlyfwd);
}

bool CarConnTwistPrim::SetPrimitives(double dt, CarPrimitiveCfg& cfg){

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
  vss.resize(grid.gs(0));
  for(size_t i=0;i<vss.size();i++)
    vss[i].reserve(cfg.na*cfg.nl);

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
            vss[idx_a].push_back(vfp);
            if(!cfg.fwdonly)
              vss[idx_a].push_back(vbp);
          }else{
            vss[idx_a].push_back(vfp);
            vss[idx_a].push_back(vfn);
            if(!cfg.fwdonly){
              vss[idx_a].push_back(vbp);
              vss[idx_a].push_back(vbn);
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
            vss[idx_a].push_back(vfp);
            if(!cfg.fwdonly)
              vss[idx_a].push_back(vbp);
          }else{
            vss[idx_a].push_back(vfp);
            vss[idx_a].push_back(vfn);
            if(!cfg.fwdonly){
              vss[idx_a].push_back(vbp);
              vss[idx_a].push_back(vbn);
            }
          }
        }
      }
    }

  }

  return true;
}

bool CarConnTwistPrim::SetPrimitives(double dt, double vx, double kmax, int kseg, bool onlyfwd) {
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

bool CarConnTwistPrim::Flow(std::tuple< SE2CellPtr, SE2Prim, double>& pathTuple,
                           const Matrix3d& g_from,
                           const Vector3d& v, bool fwd) const {

  //check if the cells encountered along the primitive and at the end, are free from obstacles or not
    double d = fabs(v[1]); // total distance along curve
    int n_seg = ceil(d/ (2 * grid.cs[1])); // 2 * grid.cs[1] is to improve efficiency
    double s = d/n_seg;
    Vector3d axy;
    Matrix3d g, dg, g_inv;
    SE2CellPtr to(nullptr);
    for (int i_seg=1; i_seg<=n_seg; i_seg++) {
      se2_exp(dg, (s*i_seg / d) * v);
      g = g_from * dg;
      se2_g2q(axy, g);
      to = grid.Get(axy);
      if (!to)
        return false;
    }
    Matrix3d g_to = to->data;

    //cost of the path
    if(fwd){
      se2_inv(g_inv,g_from);
      dg = g_inv*g_to;
    }else{
      se2_inv(g_inv,g_to);
      dg = g_inv*g_from;
    }
    Vector3d v_slip; se2_log(v_slip,dg);//twist that take you exactly to successor
    double cost = 0.1*abs(v_slip(0)) + v_slip.tail<2>().norm();//assuming carcost has ac=0.1

    //twist that takes you exactly only to (successor.x, successor.y) not angle
    Vector3d v_noslip;v_noslip << getWVx(dg(0,2),dg(1,2)),0;

    //Set the path tupule
    std::get<0>(pathTuple) = to;
    std::get<1>(pathTuple) = v_noslip;
    std::get<2>(pathTuple) = cost; //replace this by bullet prim cost if it satisfies the relationship

    return true;
}

bool CarConnTwistPrim::
operator()(const SE2Cell& from,
           std::vector< std::tuple<SE2CellPtr, SE2Prim, double> >& paths,
           bool fwd) const {
  Matrix3d g0;
  se2_q2g(g0, from.c);

  const vector<Vector3d>* p_vstemp;
  int idx_a = grid.Index(from.c,0);
  if(!vss.empty())
    p_vstemp = &(vss[idx_a]);
  else
    p_vstemp = &vs;

  paths.clear();
  //  vector< Vector3d >::const_iterator it;
  for (auto&& s : *p_vstemp) {
    // reverse time if fwd=false
    std::tuple<SE2CellPtr, SE2Prim, double> pathTuple;
    if (!Flow(pathTuple, g0, (fwd ? dt : -dt) * s,fwd))
      continue;

    assert(std::get<0>(pathTuple));


    paths.push_back(pathTuple);
  }
  return true;
}


bool CarConnTwistPrim::GetPrims(const Vector3d pos, vector<vector<Vector2d>>& prims ){
  //Display the primitive at start
  SE2CellPtr cell_start = grid.Get(pos);
  vector<std::tuple< SE2CellPtr, SE2Prim, double> > paths;
  if(cell_start){
    (*this)(*cell_start,paths,true);
    Matrix3d g0 = cell_start->data;
    prims.reserve(paths.size());
    for(auto& path:paths){
      SE2CellPtr cell_to = std::get<0>(path);
      if(!cell_to)
        continue;

      SE2Prim& v=  std::get<1>(path);
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
