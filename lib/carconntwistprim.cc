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

Vector2d xy2wt2( double xf,double yf, double u ){
  double w, t;
  if(abs(yf)<1e-12){
    w=0;
    t=xf/u;
  }else{
    w = 2*yf/(u*(xf*xf+yf*yf));
    t = atan2(w*xf/u, 1-w*yf/u)/w;
  }
  return Vector2d(w,t);
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

  cfg.nl = cfg.nl<2 ? 2: cfg.nl;  //! Make sure the number of different traj length is atleast 2
  cfg.na = cfg.na<3 ? 3: cfg.na;  //! Make sure the number of differentsteering angle is atleast 3

  double u = 1.0; //The speed can be anything. Using for clarity

  double tmin = cfg.lmin/u;
  double tmax = cfg.lmax/u;
  double wmax = u*cfg.tphioverlmax;
  double rmax = cfg.lmax/cfg.amax;
  double ymax = rmax - rmax*cos(cfg.amax);
  double xmax = cfg.lmax;

  this->dt = dt;

  double del_t = exp(log(tmax/tmin)/(cfg.nl-1));

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
      //cout<<"idx_a:"<<idx_a<<endl;
      double th = grid.xlb(0) + grid.cs(0)*(idx_a+0.5);
      for (size_t i = 0, size = grid_seen.size(); i < size; i++)
        *(grid_seen.data()+i) = false;
      Affine3d igorg_to_car = igorg_to_mgcorg * AngleAxisd(th,Vector3d::UnitZ());
      for(size_t idx_t=0 ; idx_t< cfg.nl;idx_t++){
        double t = tmin*pow(del_t, idx_t);
        for(size_t idx_w=0 ; idx_w < cfg.na; idx_w++){
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
          Vector2d wtend= xy2wt2(xyzend_snapped(0), xyzend_snapped(1),u);
          double wend = wtend(0); double tend = wtend(1);
          Vector3d vfp(wend*tend, u*tend,0);
          Vector3d vfn(-wend*tend, u*tend,0);
          Vector3d vbp(wend*tend, -u*tend,0);
          Vector3d vbn(-wend*tend, -u*tend,0);

          if(abs(wend)<1e-10){
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

    for(size_t idx_t=0 ; idx_t< cfg.nl;idx_t++){
      double t = tmin*pow(del_t, idx_t);
      for(size_t idx_w=0 ; idx_w < cfg.na; idx_w++){
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
    SE2CellPtr to(nullptr);
    Vector3d axy;
    Matrix3d g, dg,g_from_inv;
    for (int i_seg=1; i_seg<=n_seg; i_seg++) {
      se2_exp(dg, (s*i_seg / d) * v);
      g = g_from * dg;
      se2_g2q(axy, g);
      to = grid.Get(axy);
      if (!to)
        return false;
    }

    //Finding the twist element that take you exactly to the center of the cell
    se2_inv(g_from_inv,g_from);
    dg = g_from_inv*(to->data);
    Vector3d wvxvy_timesdt; se2_log(wvxvy_timesdt,dg);


    //Set the path tupule
    std::get<0>(pathTuple) = to;
    if(fwd)
      std::get<1>(pathTuple) = wvxvy_timesdt;
    else
      std::get<1>(pathTuple) = -wvxvy_timesdt;
    std::get<2>(pathTuple) = d; //replace this by bullet prim cost if it satisfies the relationship

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

}
