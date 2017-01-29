// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu> and Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_CARCONNECTIVITY_H_
#define DSL_LIB_CARCONNECTIVITY_H_

#include "gridconnectivity.h"
#include "cargrid.h"
#include "terrainse2grid.h"
#include <vector>
#include "gridcost.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace dsl {

using SE2Twist = Eigen::Vector3d;
using SE2Path = std::vector<Eigen::Matrix3d>;
using Vector1d = Eigen::Matrix<double, 1, 1>;

struct CarPrimitiveConfig {
  bool    fwdonly;      //! Decides if the car moves only in the forward direction
  double  tphioverlmax; //! Max(tan(phi)/l) possible for the car
  double  lmin;         //! Min length of the pimitive
  double  lmax;         //! Max length of the primitive
  int     nl;           //! Number of different primitive lengths from lenmin to lenmax. If 1 lmin is chosen.
  double  amax;         //! Maximum angle turned
  int     na;           //! number of steering angles from 0 to phi. If even changed to next odd number
  bool    pert;         //! preturb primitive length to have better spread of primitives
};


/**
 * Simple car connectivity using primitives. It enables generation of successors/predecessors Cells
 * given a Cell, along with the connection in between them. The connection type is a se2 twist element
 *
 * Author: Marin Kobilarov and Subhransu Mishra
 */
template<typename PointType, typename DataType, typename ConnectionType>
class SE2GridConnectivity: public GridConnectivity< PointType, DataType, ConnectionType > {
public:

  using TypedCell = Cell<PointType, DataType>;
  using TypedCellCref = typename TypedCell::Cref;
  using TypedCellPtr = typename TypedCell::Ptr;
  using TypedCellConnectionCostTuple = std::tuple<TypedCellPtr, ConnectionType, double>;
  using TypedGrid = Grid<PointType, TypedCellPtr>;
  using TypedGridCost = GridCost<PointType, DataType>;
  /**
   * Initialize cargrid connectivity
   * @param grid The car grid
   * @param cfg The configuration for generating primitives
   */
  SE2GridConnectivity(const TypedGrid& grid,const TypedGridCost& cost, CarPrimitiveConfig& cfg, bool allow_slip = true)
:grid_(grid),cost_(cost), allow_slip_(allow_slip){
    SetPrimitives(1, cfg);
  }

  /**
   * Initialize cargrid connectivity with primitives corresponding to
   * motions with constant body fixed velocities for a fixed duration dt
   * @param grid The car grid
   * @param vbs body-fixed velocities (each v=(vw,vx,vy))
   * @param dt time duration
   */
  SE2GridConnectivity(const TypedGrid& grid,const TypedGridCost& cost,
                  const std::vector<Eigen::Vector3d>& vbs,
                  double dt = .25, bool allow_slip = true) : grid_(grid),cost_(cost), vbs_(vbs), dt_(dt), allow_slip_(allow_slip)
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
  SE2GridConnectivity(const TypedGrid& grid,const TypedGridCost& cost,
                  double dt = .25,
                  double vx = 5,
                  double kmax = 0.577,
                  int kseg = 4,
                  bool onlyfwd = false,
                  bool allow_slip = true)
  : grid_(grid),cost_(cost), allow_slip_(allow_slip) {
    SetPrimitives(dt, vx, kmax, kseg, onlyfwd);
  }


  /**
   * Generates primitives based on the configuration
   * @param dt time duration
   * @param cfg the configuration for generation of primitives
   * @return
   */
  bool SetPrimitives(double dt, CarPrimitiveConfig& cfg){

    vbs_.clear();
    this->dt_ = dt;

    cfg.nl = cfg.nl<1 ? 1: cfg.nl;  //! Make sure the number of different traj length is atleast 1
    cfg.na = cfg.na<3 ? 3: cfg.na;  //! Make sure the number of differentsteering angle is atleast 3

    double u = 1.0; //The speed can be anything. Using for clarity
    double tmin = cfg.lmin/u;
    double tmax = cfg.lmax/u;
    double wmax = u*cfg.tphioverlmax;
    double del_t = exp(log(tmax/tmin)/(cfg.nl-1));
    double del_w = wmax/(cfg.na-1);

    for(int idx_t=0 ; idx_t< cfg.nl;idx_t++){ //iterate over all all lengths
      double t;
      if(cfg.nl==1)
        t = tmin;
      else
      t = tmin*pow(del_t, idx_t);
      for(int idx_w=0 ; idx_w < cfg.na; idx_w++){ //iterate over all angular velocities
        double w = idx_w*del_w;
        double tpert=t;
        if(cfg.pert && idx_t)
          tpert = t*(1 + 0.2*(cos(40*w)-1)); // t perturbed
        else
          tpert = t;
        double apert = w*tpert;              // angle perturbed
        double lpert = u*tpert;
        if(apert > cfg.amax || lpert > 1.5*cfg.lmax || lpert < 0.6*cfg.lmin)
          continue;

        Eigen::Vector3d vfp(w*tpert, u*tpert,0);  //forward vel positive angular vel
        Eigen::Vector3d vfn(-w*tpert, u*tpert,0); //forward vel negative angular vel
        Eigen::Vector3d vbp(w*tpert, -u*tpert,0); //reverse vel positive angular vel
        Eigen::Vector3d vbn(-w*tpert, -u*tpert,0);//reverse vel negative angular vel

        if(std::fabs(w)<1e-16){
          vbs_.push_back(vfp);
          if(!cfg.fwdonly)
            vbs_.push_back(vbp);
        }else{
          vbs_.push_back(vfp);
          vbs_.push_back(vfn);
          if(!cfg.fwdonly){
            vbs_.push_back(vbp);
            vbs_.push_back(vbn);
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

    this->dt_ = dt;

    vbs_.clear();

    for (int i = 0; i <= kseg; i++) {
      double k = i*kmax/kseg;
      double w = vx*k;
      vbs_.push_back(Eigen::Vector3d(w, vx, 0));
      vbs_.push_back(Eigen::Vector3d(-w, vx, 0));
      if (!onlyfwd) {
        vbs_.push_back(Eigen::Vector3d(w, -vx, 0));
        vbs_.push_back(Eigen::Vector3d(-w, -vx, 0));
      }
    }

    return true;
  }

  bool operator()(const TypedCell& from,
                  std::vector<TypedCellConnectionCostTuple>& paths,
                  bool fwd = true) const override;

  bool Free(const DataType &g) const override {
    return true;
  }

  /**
   * Utility function to get all the primitives starting at a given position. Each primitive is a vector of
   * positions along the path. This is particularly useful for plotting the primitives.
   * @param prims A single point along a primitive is xy pos
   */
  bool GetPrims(const Eigen::Vector3d& pos, std::vector<std::vector<Eigen::Vector2d>>& prims ){
    //Display the primitive at start
    TypedCellPtr cell_start = grid_.Get(pos);
    if(!cell_start){
      prims.clear();
      return false;
    }
    std::vector<TypedCellConnectionCostTuple> paths;
    if(cell_start){
      (*this)(*cell_start,paths,true);

      Eigen::Matrix3d g0;
      se2_q2g(g0, cell_start->c);

      prims.reserve(paths.size());
      for(auto& path:paths){
        TypedCellPtr cell_to = std::get<0>(path);
        if(!cell_to)
          continue;
        Eigen::Matrix3d gto;
        se2_q2g(gto, cell_to->c);

        Eigen::Matrix3d gi,dg;
        se2_inv(gi,g0);
        dg = gi*gto; //relative of from.data to to->data

        Eigen::Vector3d v;
        se2_log(v,dg);//twist that take you exactly to successor

        double d = std::fabs(v[1]); // total distance along curve
        int n_seg = ceil(d/ (2 * grid_.cs()[1])); // 2 * grid.cs[1] is to improve efficiency
        double s = d/n_seg;
        std::vector<Eigen::Vector2d> prim(0);
        prim.reserve(n_seg+1);
        prim.push_back(Eigen::Vector2d(g0(0,2),g0(1,2)));
        for (int i_seg=1; i_seg<=n_seg; i_seg++) {
          Eigen::Vector3d axy;
          Eigen::Matrix3d g, dg;
          se2_exp(dg, (s*i_seg / d) * v);
          g = g0 * dg;
          prim.push_back(Eigen::Vector2d(g(0,2),g(1,2)));
        }
        prims.push_back(prim);
      }
      return true;
    }else{
      return false;
    }
  }

  const TypedGrid& grid_; ///< the grid
  const TypedGridCost& cost_; ///< cost interface that gives you cost of taking primitive and if path is blocked
  std::vector< Eigen::Vector3d > vbs_; ///< primitives defined using motions with constant body-fixed velocities (w,vx,vy)
  double dt_ = .5; ///< how long are the primitives
  bool allow_slip_; ///< to use primitive with slip or the no_slip counterpart? slip is accurate but not suitable for car
};

using CarTwistConnectivity     = SE2GridConnectivity<SE2Cell::PointType, SE2Cell::DataType, SE2Twist>;
using CarPathConnectivity      = SE2GridConnectivity<SE2Cell::PointType, SE2Cell::DataType, SE2Path>;
using TerrainTwistConnectivity = SE2GridConnectivity<SE2TerrainCell::PointType, SE2TerrainCell::DataType, SE2Twist>;
}

#endif //DSL_CARCONNECTIVITY_H
