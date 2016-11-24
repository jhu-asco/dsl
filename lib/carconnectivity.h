// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu> and Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARCONNECTIVITY_H
#define DSL_CARCONNECTIVITY_H

#include "carprimitiveconfig.h"
#include "gridconnectivity.h"
#include "cargrid.h"
#include <vector>
#include "carcost.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace dsl {

using SE2Twist = Eigen::Vector3d;
using SE2Path = std::vector<Eigen::Matrix3d>;
using Vector1d = Eigen::Matrix<double, 1, 1>;

/**
 * Simple car connectivity using primitives. It enables generation of successors/predecessors Cells
 * given a Cell, along with the connection in between them. The connection type is a se2 twist element
 *
 * Author: Marin Kobilarov and Subhransu Mishra
 */
template<typename ConnectionType>
class CarConnectivity: public GridConnectivity< SE2Cell::PointType, SE2Cell::DataType, ConnectionType > {
public:

  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param cfg The configuration for generating primitives
   */
  CarConnectivity(const CarGrid& grid,const CarCost& cost, CarPrimitiveConfig& cfg)
:grid(grid),cost(cost){
    SetPrimitives(1, cfg);
  }

  /**
   * Initialize cargrid connectivity with primitives corresponding to
   * motions with constant body fixed velocities for a fixed duration dt
   * @param grid The car grid
   * @param vbs body-fixed velocities (each v=(vw,vx,vy))
   * @param dt time duration
   */
  CarConnectivity(const CarGrid& grid,const CarCost& cost,
                  const std::vector<Eigen::Vector3d>& vbs,
                  double dt = .25) : grid(grid),cost(cost), vbs(vbs), dt(dt)
  {
  }

  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param dt time duration
   * @param vx forward velocity.
   * @param kmax maximum curvature k = Tan(max steering angle)/axle_length; Default value is 0.577=tan(M_PI/6)/1
   * @param kseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param onlyfwd If true, then only +ve forward velocity is used
   */
  CarConnectivity(const CarGrid& grid,const CarCost& cost,
                  double dt = .25,
                  double vx = 5,
                  double kmax = 0.577,
                  int kseg = 4,
                  bool onlyfwd = false)
  : grid(grid),cost(cost) {
    SetPrimitives(dt, vx, kmax, kseg, onlyfwd);
  }


  /**
   * Generates primitives based on the configuration
   * @param dt time duration
   * @param cfg the configuration for generation of primitives
   * @return
   */
  bool SetPrimitives(double dt, CarPrimitiveConfig& cfg){

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

  /**
   * Use a set of primitive motions, i.e. arcs with body fixed forward
   * velocities (-v,v) and
   * angular velocity (-w,0,w), lasting time duration dt. By default there are 6
   * such combinations.
   * In addition we made it configurable so that more primitves with angular
   * velocities in between(-w,w)
   * can be added in between by choosing the kseg parameter >2. Also only +ve
   * forward velocity can be chosen
   * @param dt time duration
   * @param vx forward velocity.
   * @param kmax maximum curvature k = Tan(max steering angle)/axle_length; Default value is 0.577=tan(M_PI/6)/1
   * @param kseg It decides the discretization of max angular velocity when
   * making the motion primitives
   * @param onlyfwd If true, then only +ve forward velocity is used
   * @return true on success
   */
  bool SetPrimitives(double dt,
                     double vx,
                     double kmax,
                     int kseg = 4,
                     bool onlyfwd = false){
    if (dt <= 0)
      return false;

    kseg = kseg < 1 ? 1 : kseg;

    this->dt = dt;

    vbs.clear();

    for (int i = 0; i <= kseg; i++) {
      double k = i*kmax/kseg;
      double w = vx*k;
      vbs.push_back(Vector3d(w, vx, 0));
      vbs.push_back(Vector3d(-w, vx, 0));
      if (!onlyfwd) {
        vbs.push_back(Vector3d(w, -vx, 0));
        vbs.push_back(Vector3d(-w, -vx, 0));
      }
    }

    return true;
  }

  bool operator()(const SE2Cell& from,
                  std::vector< std::tuple<SE2CellPtr, ConnectionType, double> >& paths,
                  bool fwd = true) const override;

  bool Free(const Eigen::Matrix3d &g) const override {
    return true;
  }

  /**
   * Utility function to get all the primitives starting at pos. Each primitive is a vector of
   * positions along the path
   * @param prims A single point along a primitive is xy pos
   */
  bool GetPrims(const Vector3d pos, vector<vector<Vector2d>>& prims ){
    //Display the primitive at start
    SE2CellPtr cell_start = grid.Get(pos);
    if(!cell_start){
      prims.clear();
      return false;
    }
    vector<std::tuple< SE2CellPtr, ConnectionType, double> > paths;
    if(cell_start){
      (*this)(*cell_start,paths,true);
      Matrix3d g0 = cell_start->data;
      prims.reserve(paths.size());
      for(auto& path:paths){
        SE2CellPtr cell_to = std::get<0>(path);
        if(!cell_to)
          continue;

        Matrix3d gi,dg; se2_inv(gi,g0); dg = gi*cell_to->data; //relative of from.data to to->data
        Vector3d v; se2_log(v,dg);//twist that take you exactly to successor
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


  const CarGrid& grid; ///< the grid

  std::vector< Eigen::Vector3d > vbs; ///< primitives defined using motions with constant body-fixed velocities (w,vx,vy)
  std::vector< std::vector<Eigen::Vector3d > > vbs_per_angle; ///< A set of primitives for each angle

  double dt = .5; ///< how long are the primitives

  const CarCost& cost;
};

using CarTwistConnectivity = CarConnectivity<SE2Twist>;
using CarPathConnectivity  = CarConnectivity<SE2Path>;
}

#endif //DSL_CARCONNECTIVITY_H
