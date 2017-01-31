// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu> and Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_GRID_H_
#define DSL_LIB_GRID_H_

#include "cell.h"
#include "utils.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <exception>
#include <type_traits>
#include "grid_data.pb.h"
#include "google/protobuf/io/coded_stream.h"

namespace dsl {

struct EmptyData {};
/**
 * An n-dimenensional grid consisting of abstract "cells", or elements
 * identified by a set of coordinates of type PointType. The cells can hold
 * a specified DataType or a unique_ptr to the DataType. The second option
 * allows for allocation of data for specific cells. All that is controlled by
 * the template parameter UsePtr
 *
 * Explanation of certain important variables:
 *   id:   id for direct lookup in an array(or std::vector)
 *   i:    dimension number(for grid over xyz coordinates i=0 refers to the x dimension)
 *   gidx: grid index to locate a cell in a grid. For a 3d grid, gidx=(0,0,0) is the first cell along all dimension.
 *   idx:  grid index along a particular dimension
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 *
 *\todo: 1. Make sure origin lies at center of a cell
 *       2. Check for bound in GetSlice and GetStack methods
 *       3. Implement GetStack method that can perform occupancy dilation
 *       4. Add clone methods to deep copy data
 *
 */
template < class PointType, class DataType, bool UsePtr = true>
class Grid {
public:

  //n-dimensional std::vectors
  using Vectornd   =  Eigen::Matrix< double,  PointType::SizeAtCompileTime,   1 >;
  using Vectorni   =  Eigen::Matrix< int,     PointType::SizeAtCompileTime,   1 >;
  using Vectornb   =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime,   1 >;

  //n-1 dimensional std::vectors
  using Vectornm1d =  Eigen::Matrix< double,  PointType::SizeAtCompileTime-1, 1 >;
  using Vectornm1i =  Eigen::Matrix< int,     PointType::SizeAtCompileTime-1, 1 >;
  using Vectornm1b =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime-1, 1 >;

  //n+1 dimensional std::vectors
  using Vectornp1d =  Eigen::Matrix< double,  PointType::SizeAtCompileTime+1, 1 >;
  using Vectornp1i =  Eigen::Matrix< int,     PointType::SizeAtCompileTime+1, 1 >;
  using Vectornp1b =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime+1, 1 >;

  using Ptr = std::shared_ptr<Grid>;

  using Slice = Grid<Vectornm1d, DataType, UsePtr>;
  using SlicePtr = std::shared_ptr<Slice>;

  using Stack = Grid<Vectornp1d, DataType, UsePtr>;
  using StackPtr = std::shared_ptr<Stack>;

  using CellType = typename std::conditional<UsePtr, std::unique_ptr<DataType>, DataType>::type;
  static const std::integral_constant<bool, UsePtr> cells_store_ptr_type; //std::false_type or true_type
  static const bool cells_store_ptr = UsePtr;

  /**
   * Initialize the grid using state lower bound, state upper bound, the number of grid cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of grid cells per each dimension
   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   */
  Grid(const Vectornd& xlb, const Vectornd& xub, const Vectorni& gs, const Vectornb& wd = Vectornb::Zero())
  : xlb_(xlb), xub_(xub), gs_(gs), wd_(wd), cells_(0){
    ds_ = xub - xlb;
    nc_ = 1;
    cgs_[0] = 1;
    for (int i = 0; i < n_; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      nc_ *= gs[i]; // total number of cells
      cs_[i] = ds_[i] / gs[i];
      if(i>0)
        cgs_[i] = cgs_[i-1]*gs[i-1];
    }

    cells_.resize(nc_);
  }

  //  /**
  //   * Initialize the map using state lower bound, state upper bound, and a suggested cell size
  //   * @param xlb state lower bound
  //   * @param xub state upper bound
  //   * @param scs suggested cell size along each dimension. The true cs is the one which results in
  //   * integer number of grid cells between xlb and xub
  //   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
  //   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
  //   */
  //  Grid(const Vectornd& xlb, const Vectornd& xub, const Vectornd& scs, const Vectornb wd = Vectornb::Zero())
  //  : n(xlb.size()), xlb(xlb), xub(xub), wd(wd), cells(0){
  //    ds = xub - xlb;
  //    nc = 1;
  //    cgs[0] = 1;
  //    for (int i = 0; i < n; ++i) {
  //      assert(xlb[i] < xub[i]);
  //      assert(scs[i] > 0);
  //      gs[i] = floor(ds[i] / scs[i]);
  //      cs[i] = ds[i] / gs[i];
  //      nc *= gs[i]; // total number of cells
  //      if(i>0)
  //        cgs[i] = cgs[i-1]*gs[i-1];
  //    }
  //    cells.resize(nc);
  //  }

  /**
   * Initialize the map using a suggested state lower bound, suggested state upper bound, and a suggested cell size.
   *
   * @param sxlb state lower bound
   * @param sxub state upper bound
   * @param scs suggested cell size along each dimension. The true cs is the one which results in
   * integer number of grid cells between xlb and xub
   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   * @param ss Special Settings for each dimension. If set to true for a wrapped dimension then it the number of cell are
   * multiples of 4 and the xlb and xub are shifted by half cell size. This is useful for making sure, in case of angles,
   * that when moving forward angles don't change if starting from cells with yaw = -pi,-pi/2, 0, pi/2. If set to true for
   * flat dimensions it makes sure that 0 lies at a cell center of a cell in the grid( or extrapolated grid if center is not
   * contained in the grid)
   */
  Grid(const Vectornd& sxlb, const Vectornd& sxub, const Vectornd& scs,
       const Vectornb& wd = Vectornb::Zero(), const Vectornb& ss = Vectornb::Zero())
  : wd_(wd), cells_(0){
    nc_ = 1;
    cgs_[0] = 1;
    for (int i = 0; i < n_; ++i) {
      assert(sxlb[i] < sxub[i]);
      assert(scs[i] > 0);

      if(ss[i]){
        if(!wd[i]){ //if dimension is flat(cs=scs)
          cs_[i] = scs[i];
          int n_org2ub = floor(sxub[i]/cs_[i] - 0.5);
          int n_org2lb =  ceil(sxlb[i]/cs_[i] + 0.5);
          gs_[i] = n_org2ub  - n_org2lb + 1;
          xub_[i] = (n_org2ub +0.5)*cs_[i];
          xlb_[i] = (n_org2lb -0.5)*cs_[i];
          ds_[i] = xub_[i] - xlb_[i];
        }else{ //dimension is wrapped
          ds_[i] = sxub[i] - sxlb[i]; //For angles this is = 2*M_PI
          gs_[i] = 4*ceil(ds_[i]/(4*scs[i])); //Multiple of 4 number of cells
          cs_[i] = ds_[i]/gs_[i];
          xlb_[i] = -ds_[i]/2 + cs_[i]/2; //shifted by half cell so that pi/(pi/2)/0 lies at center of a cell
          xub_[i] =  ds_[i]/2 + cs_[i]/2; //shifted by half cell so that pi/(pi/2)/0 lies at center of a cell
        }
      }else{
        xlb_[i] = sxlb[i];
        xub_[i] = sxub[i];
        ds_[i] = xub_[i] - xlb_[i];
        gs_[i] = floor(ds_[i] / scs[i]);
        cs_[i] = ds_[i] / gs_[i];
      }

      nc_ *= gs_[i]; // total number of cells
      if(i>0)
        cgs_[i] = cgs_[i-1]*gs_[i-1];
    }
    cells_.resize(nc_);
  }

  /**
   * Copies the grid structure and allocates memory for all but assigns only if cell doesn't store shared_ptr
   * @param gridcore
   */
  Grid(const Grid &gridcore)
  : nc_(gridcore.nc_), xlb_(gridcore.xlb_), xub_(gridcore.xub_), ds_(gridcore.ds_),
    cs_(gridcore.cs_), gs_(gridcore.gs_), cgs_(gridcore.cgs_), wd_(gridcore.wd_){
    cells_.resize(nc_);

    //copy data if not shared_ptr
    if(!cells_store_ptr)
      cells_ = gridcore.cells_;
  }

  virtual ~Grid() {}

  /**
   * Copies the grid structure and allocates memory for all but assigns only for arithmetic CellType
   * @param The grid to copy the structure from
   * @return
   */
  Grid& operator=(const Grid& gridcore){
    nc_  = gridcore.nc_;
    xlb_ = gridcore.xlb_;
    xub_ = gridcore.xub_;
    gs_  = gridcore.gs_;
    cs_  = gridcore.cs_;
    cgs_ = gridcore.cgs_;
    wd_  = gridcore.wd_;
    cells_.resize(nc_);

    //copy data if not shared_ptr
    if(!cells_store_ptr)
      cells_ = gridcore.cells_;
    return *this;
  }

  /**
   * Stacks up the grid along a new dimension. Say the current grid
   * has 3 dimensions(angle/a, x and y), a new dimension is inserted as:
   *    a   x   y
   *  ^   ^   ^   ^
   *  0   1   2   3
   * @param i coordinated index/ dimension
   * @param lbi lower bound along that dimension
   * @param ubi upper bound along that dimension
   * @param gsi grid size along that dim
   * @param wdi dimension is wrapped or not
   * @param init initialize the Stack from grid
   * @return
   */
  StackPtr GetStack(int i, double lbi, double ubi, int gsi, bool wdi = false, bool init = true) const {
    Vectornp1d xlb_stack, xub_stack;
    Vectornp1i gs_stack;
    Vectornp1b wd_stack;

    //Get the new dimension
    xlb_stack = InsertDimension(xlb_,i, lbi);
    xub_stack = InsertDimension(xub_,i, ubi);
    gs_stack  = InsertDimension(gs_, i, gsi);
    wd_stack  = InsertDimension(wd_, i, wdi);

    StackPtr pstack(new Stack(xlb_stack,xub_stack, gs_stack, wd_stack));

    if(!init || cells_store_ptr)
      return pstack;

    //Iterate over all cells
    auto fun = [&](int id, const Vectorni& gidx) {
      for( int idx=0; idx < pstack->gs_[i] ; idx++ ){ //repeat val it along the extra dimension
        Vectornp1i gidx_stack = InsertDimension(gidx,i,idx);
        pstack->cells_[pstack->Id(gidx_stack)] = cells_[id]; //get value of that cell
      }
    };
    LoopOver(fun);

    return pstack;
  }

  /**
   * Stacks up the grid along a new dimension. Say the current grid
   * has 3 dimensions(a, x and y), a new dimension is inserted as:
   *    a   x   y
   *  ^   ^   ^   ^
   *  0   1   2   3
   * @param i coordinated index/ dimension
   * @param lbi lower bound along that dimension
   * @param ubi upper bound along that dimension
   * @param scsi suggested cell size along that dim
   * @param wdi dimension is wrapped or not
   * @param init initialize the Stack from grid
   * @return
   */
  StackPtr GetStack(int i, double lbi, double ubi, double scsi, bool wdi = false, bool init = true) const {
    Vectornp1d xlb_stack, xub_stack,scs_stack;
    Vectornp1b wd_stack;


    //Get the new dimension
    xlb_stack = InsertDimension(xlb_,i, lbi);
    xub_stack = InsertDimension(xub_,i, ubi);
    scs_stack = InsertDimension(cs_, i, scsi);
    wd_stack  = InsertDimension(wd_, i, wdi);

    StackPtr pstack(new Stack(xlb_stack,xub_stack, scs_stack, wd_stack));

    //copy data if not shared_ptr
    if(!init || cells_store_ptr)
      return pstack;

    for(int id=0; id < nc_; id++){
      Vectorni gidx;
      Index(id,&gidx);
      CellType val = cells_[id];
      for( int idx=0; idx < pstack->gs()[i] ; idx++ ){
        Vectornp1i gidx_stack = InsertDimension(gidx,i,idx);
        pstack->cells_[pstack->Id(gidx_stack)] = val;
      }
    }
    return pstack;
  }

  /**
   * Create a slice of the current grid along dimension at particular index
   * @param idx The dimension along which the slice is created
   * @param i coordinated index/ dimension
   * @param init initialize the slice from grid
   * @return a shared_ptr to the grid slice created
   */
  SlicePtr GetSlice(int idx, int i, bool init = true) const {
    Vectornm1d xlb_slice, xub_slice;
    Vectornm1i gs_slice;
    Vectornm1b wd_slice;

    //Get the new dimension
    xlb_slice = RemoveDimension(xlb_,i);
    xub_slice = RemoveDimension(xub_,i);
    gs_slice  = RemoveDimension(gs_,i);
    wd_slice  = RemoveDimension(wd_,i);

    SlicePtr pslice(new Slice(xlb_slice,xub_slice, gs_slice, wd_slice));

    if(!init || cells_store_ptr)
      return pslice;

    for(int id_slice=0; id_slice < pslice->nc_; id_slice++){
      Vectornm1i gidx_slice;
      pslice->Index(id_slice, &gidx_slice);
      Vectorni gidx = InsertDimension(gidx_slice,i,idx);
      pslice->cells_[id_slice] = Get(Id(gidx));
    }
    return pslice;
  }

  /**
   * Create a slice of the current grid along dimension at particular index
   * @param val val = point[i] where i is the dimension along which the slice is created
   * @param i coordinated index/ dimension
   * @param init initialize the slice from grid
   * @return a shared_ptr to the grid slice created
   */
  SlicePtr GetSlice(double val, int i, bool init = true) const {
    int idx = Index(val,i);
    return GetSlice(idx,i,init);
  }

  /**
   * Create a slice of the current grid along dimension at particular index
   * @param pslice pointer to the slice
   * @param idx index in a grid along the dim coordinate/dimension
   * @param dim coordinated index/ dimension
   */
  void GetSlice(Slice* pslice, int idx, int dim, bool init = true) const {
    Slice& slice = *pslice;
    //Get the new dimension
    slice.xlb_ = RemoveDimension(xlb_,dim);
    slice.xub_ = RemoveDimension(xub_,dim);
    slice.gs_  = RemoveDimension(gs_,dim);
    slice.wd_  = RemoveDimension(wd_,dim);
    slice.cs_  = RemoveDimension(cs_,dim);
    slice.ds_  = RemoveDimension(ds_,dim);
    slice.nc_=1;
    slice.cgs_[0]=1;
    for(int i = 0; i < n_- 1; i++){
      slice.nc_ *=slice.gs_[i];
      if(i>0)
        slice.cgs_[i] = slice.cgs_[i-1]*slice.gs_[i-1];
    }
    slice.cells_.resize(slice.nc());//if data is already allocated resize does nothing

    if(!init || cells_store_ptr)
      return;

    for(int id_slice=0; id_slice < slice.nc(); id_slice++){
      Vectornm1i midx_slice;
      slice.Index(id_slice, &midx_slice);
      Vectorni midx = InsertDimension(midx_slice,dim,idx); //midx is multidim index
      slice.cells_[id_slice] = Get(Id(midx));
    }
  }


  /**
   * Create new Grid object whose cell resolution is scale times higher.
   * Useful for increasing resolution of a map for display
   * @param scale
   * @param init initialize the scaled map from grid
   * @return
   * TODO: Add option of linear/quadratic/cubic interpolation when scaling
   */
  Ptr ScaleUp(int scale, bool init = true) const {
    Ptr pscaled;
    if(scale<1)
      return pscaled;

    if(scale==1){
      pscaled.reset(new Grid(*this));
      return pscaled;
    }

    Vectorni gs_scaled = gs_*scale;
    pscaled.reset(new Grid(xlb_,xub_, gs_scaled, wd_));

    if(!init || cells_store_ptr)
      return pscaled;

    for(int id_scaled = 0; id_scaled <pscaled->nc(); id_scaled++){
      Vectornd cc;
      pscaled->CellCenter(id_scaled, &cc);
      pscaled->cells_[id_scaled] = Get(cc,false);
    }
    return pscaled;
  }


  //  /**
  //   * Takes a set of vertices which are corners of a kernel map. The kernel.cs == grid.cs by design
  //   * @param grid the main grid
  //   * @param vertices set of vertices
  //   * @param valin value inside the area marked by the vertices
  //   * @param valout value outside the area marked by the vertices
  //   */
  //  Grid::Ptr GetKernel(const std::vector<Vectornd>& vertices,CellType valin, CellType valout) const{
  //
  //    //convert vertices to grid coordinates
  //    std::vector<Vectornd> verts_grid; ToGridCoordinates(verts_grid,vertices);
  //
  //    //Get center of points and convert to grid coords
  //    Vectornd center = std::accumulate(vertices.begin(),vertices.end(),Vectornd::Zero().eval())/vertices.size();
  //    Vectornd center_grid; ToGridCoordinates(center_grid,center);
  //
  //    //Snap the vertices to cell centers while making sure it is expanded outwards
  //    std::vector<Vectorni> verts_snapped;
  //    for(size_t i=0; i<vertices.size(); i++){
  //      Vectornd verts_grid_centered = verts_grid[i] - center_grid;
  //      for(int dim=0; dim < n; dim++){
  //        verts_snapped[i][dim] = verts_grid_centered[dim]>0 ?
  //                                ceil(verts_grid[i][dim] ):floor(verts_grid[i][dim] );
  //      }
  //    }
  //
  //
  //
  //    Grid::Ptr kernel;
  //    return kernel;
  //  }

  /**
   * Check if point x is within grid bounds
   * @param x point
   * @param eps min distance from boundary for validity(only for flat dimensions)
   * @return true if within bounds
   */
  bool Valid(const Vectornd& x, double eps = 1e-10) const {
    for (int i = 0; i < x.size(); ++i) {
      if(!wd_[i]){ //if dimension is flat
        if (x[i] < xlb_[i] + eps)
          return false;
        if (x[i] > xub_[i] - eps)
          return false;
      }
    }
    return true;
  }

  /**
   * Check if a multidimensional index is within grid bounds
   * @param gidx grid index
   * @return true if within bounds
   */
  bool Valid(const Vectorni& gidx) const{
    for(size_t i=0; i< n_; i++){
      if(gidx[i]<0 || gidx[i]>=gs_[i])
        return false;
    }
    return true;
  }

  /**
   * Check if id( array lookup id) is valid
   * @param id fast lookup id on array
   * @return true if within bounds
   */
  bool Valid(int id) const{
    return (id>=0 && id<nc_);
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param gidx Multidimensional index (can extend outside the grid)
   * @param x The center point (meaningful even when gidx outside bounds)
   * @return true if cell center is inside grid bounds
   */
  bool CellCenter(const Vectorni& gidx, Vectornd* x) const{
    *x = xlb_.array() +  (gidx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs_.array() ;
    return Valid(gidx);
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param gidx Multidimensional index (can extend outside the grid)
   * @return The center point (meaningful even when gidx outside bounds)
   */
  Vectornd CellCenter(const Vectorni& gidx) const{
    Vectornd x = xlb_.array() +  (gidx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs_.array() ;
    return x;
  }

  /**
   * Finds the CellCenter corresponding to a point x_in
   * @param x_in a point for which a snapped-in cell center is needed
   * @param x The center point (meaningful even when gidx outside bounds)
   * @return true if cell center inside grid bounds
   */
  bool CellCenter(const Vectornd& x_in, Vectornd* x) const{
    Vectorni gidx;
    Index(x_in, &gidx);
    *x = xlb_.array() +  (gidx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs_.array() ;
    return Valid(gidx);
  }

  /**
   * Finds the CellCenter corresponding to a point x_in
   * @param x_in a point for which a snapped-in cell center is needed
   * @return The center point (meaningful even when gidx outside bounds)
   */
  Vectornd CellCenter(const Vectornd& x_in) const{
    Vectorni gidx;
    Index(x_in, &gidx);
    Vectornd x = xlb_.array() +  (gidx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs_.array() ;
    return x;
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param id fast access id in the cell array
   * @param x The updated center ( set to a std::vector of NANs when id not valid)
   * @return true if id is valid else false
   */
  bool CellCenter(int id, Vectornd* x) const{
    Vectorni gidx;
    if(!Index(id, &gidx)){
      *x = Vectornd::Constant(std::numeric_limits<double>::quiet_NaN());
      return false;
    }
    *x = xlb_.array() +  (gidx.template cast<double>()+Vectornd::Constant(0.5)).array()*cs_.array() ;
    return true;
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param id fast access id in the cell array
   * @return x The updated center ( set to a std::vector of NANs when id not valid)
   */
  Vectornd CellCenter(int id) const{
    Vectornd x;
    Vectorni gidx;
    if(!Index(gidx,id)){
      x = Vectornd::Constant(std::numeric_limits<double>::quiet_NaN());
      return false;
    }
    x = xlb_.array() +  (gidx.template cast<double>()+Vectornd::Constant(0.5)).array()*cs_.array() ;
    return x;
  }

  /**
   * Returns cc[i] where cc is the cell center of any cell with gidx[i] = idx
   * @param idx gidx[i] = idx, i.e. grid index along dimension i.
   * @param i coordinated index/ dimension
   * @return Returns cc[i]. Meaningful result even when idx out of range
   */
  double CellCenterIth(int idx, int i) const{
    return xlb_[i] + (idx +0.5)*cs_[i];
  }

  /**
   * Get id of point x useful for direct lookup in the cell array
   * @param x point
   * @return the cell array id. -1 if point x is not in grid
   */
  int Id(const Vectornd& x) const {
//    Vectorni gidx;
//    Index(x, &gidx);
//
//    if(Valid(gidx))
//      return gidx.transpose()*cgs_;
//    else
//      return -1;

    //This is forloop is actually faster than unrolled because there
    //are no function calls.
    int id = 0;
    for (int i = 0; i < n_; ++i) {
      // index of i-th dimension
      double xi = x[i];
      if(xi - xlb_[i] >= 0)
        xi =  fmod(xi - xlb_[i], ds_[i]) + xlb_[i];
      else
        xi = -fmod(xlb_[i] - xi, ds_[i]) - xlb_[i];
      int ind = floor((xi - xlb_[i]) / ds_[i] * gs_[i]);
      if(ind < 0 || ind >= gs_[i])
        return -1;
      id += cgs_[i] * ind;
    }
    return id;
  }

  /**
   * Get the id of a point corresponding to the grid index
   * @param gidx grid index
   * @return the cell array id. -1 if gidx is not within bounds
   */
  int Id(const Vectorni& gidx) const {
    if(Valid(gidx))
      return gidx.dot(cgs_);
    else
      return -1;
  }

  /**
   * Get the grid index of the point x along i coordinate index.
   * @param x point
   * @param i coordinated index/ dimension
   * @return gidx[i], i.e. grid index along i coordinate/dimension. Meaningful even when x is out of bounds.
   */
  int Index(const Vectornd& x, int i) const {
    double xi = x[i];
    if(wd_[i]){ //dimension is wrapped
      if(xi - xlb_[i] >= 0)
        xi =  fmod(xi - xlb_[i], ds_[i]) + xlb_[i];
      else
        xi = -fmod(xlb_[i] - xi, ds_[i]) - xlb_[i];
    }
    return floor((xi - xlb_[i]) / ds_[i] * gs_[i]);
  }

  /**
   * Get the grid index of i coordinate/dimension
   * @param xi value of a point at i coordinate/dimension
   * @param i coordinated index/ dimension
   * @return gidx[i], i.e. grid index along i coordinate/dimension. Meaningful even when x is out of bounds.
   */
  int Index(double xi, int i) const {
    if(wd_[i]){ //dimension is wrapped
      if(xi - xlb_[i] >= 0)
        xi =  fmod(xi - xlb_[i], ds_[i]) + xlb_[i];
      else
        xi = -fmod(xlb_[i] - xi, ds_[i]) - xlb_[i];
    }
    return floor((xi - xlb_[i]) / ds_[i] * gs_[i]);
  }

  /**
   * Get the grid index along all the dimensions
   * @param x point
   * @param gidx grid index. Meaningful even if x is not inside grid
   */
  void Index(const Vectornd& x, Vectorni* gidx) const {
    for(size_t i = 0; i < n_; i++){
      double xi = x[i];
      if(wd_[i]){ //dimension is wrapped
        if(xi - xlb_[i] >= 0)
          xi =  fmod(xi - xlb_[i], ds_[i]) + xlb_[i];
        else
          xi = -fmod(xlb_[i] - xi, ds_[i]) - xlb_[i];
      }
      (*gidx)[i] = floor((xi - xlb_[i]) / ds_[i] * gs_[i]);
    }
  }

  /**
   * Get the grid index along all the dimensions, for a given point
   * @param x point
   * @return grid index. Meaningful even if x is not inside grid
   */
  Vectorni Index(const Vectornd& x) const {
    Vectorni& gidx;
    for(size_t i=0;i<n_;i++){
      double xi = x[i];
      if(wd_[i]){ //dimension is wrapped
        if(xi - xlb_[i] >= 0)
          xi =  fmod(xi - xlb_[i], ds_[i]) + xlb_[i];
        else
          xi = -fmod(xlb_[i] - xi, ds_[i]) - xlb_[i];
      }
      gidx[i] = floor((xi - xlb_[i]) / ds_[i] * gs_[i]);
    }
    return gidx;
  }

  /**
   * Get the grid index along all the dimensions from the direct lookup id
   * @param id id for direct lookup in the grid array
   * @param gidx updated grid index. Set to NANs if id out of range
   * @return false if index is out of range
   */
  bool Index(int id, Vectorni* gidx) const {
    Vectorni& gidx2 = *gidx;
    if(id>=nc_ || id<0){
      gidx2 = Vectorni::Constant(std::numeric_limits<double>::quiet_NaN());
      return false;
    }
    switch(n_){
      case 1:
        gidx2[0] = id;
        break;
      case 2:
        gidx2[1] = int(id/cgs_[1]);
        gidx2[0] = int(id - gidx2[1]*cgs_[1]);
        break;
      case 3:
        gidx2[2] = int(id/cgs_[2]);
        gidx2[1] = int(id/cgs_[1] - gidx2[2]*cgs_[2]/cgs_[1]);
        gidx2[0] = int(id - gidx2[2]*cgs_[2] - gidx2[1]*cgs_[1]);
        break;
      case 4:
        gidx2[3] = int( id/cgs_[3]);
        gidx2[2] = int((id - gidx2[3]*cgs_[3])/cgs_[2]);
        gidx2[1] = int((id - gidx2[3]*cgs_[3] - gidx2[2]*cgs_[2])/cgs_[1]);
        gidx2[0] = int( id - gidx2[3]*cgs_[3] - gidx2[2]*cgs_[2] - gidx2[1]*cgs_[1]);
        break;
      case 5:
        gidx2[4] = int( id/cgs_[4]);
        gidx2[3] = int((id - gidx2[4]*cgs_[4])/cgs_[3]);
        gidx2[2] = int((id - gidx2[4]*cgs_[4] - gidx2[3]*cgs_[3])/cgs_[2]);
        gidx2[1] = int((id - gidx2[4]*cgs_[4] - gidx2[3]*cgs_[3] - gidx2[2]*cgs_[2])/cgs_[1]);
        gidx2[0] = int( id - gidx2[4]*cgs_[4] - gidx2[3]*cgs_[3] - gidx2[2]*cgs_[2] - gidx2[1]*cgs_[1]);
        break;
      default:
        for(int i=n_-1; i>=0; i--){
          gidx2[i] = int(id/cgs_[i]);
          id -= gidx2[i]*cgs_[i];
        }
    }
    return true;
  }

  /**
   * Get const ref to *cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return const ref to *cell at position x
   * If cell is unallocated then a null reference is returned
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<UsesPtr>::type >
  const DataType* Get(const Vectornd& x, bool checkValid = true) const {
    if (checkValid)
      if (!Valid(x))
        return nullptr;

    return cells_[Id(x)].get();
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return const ref to contents of cell.
   * If cell holds pointer and is not allocate then a null reference is returned
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<UsesPtr>::type >
  const DataType* Get(int id) const {
    if (id<0 || id >= nc_)
      return nullptr;

    return cells_[id].get();
  }

  /**
   * Get the cell at a given grid index
   * @param gidx grid index
   * @return Copy of contents of cell, could be shared_ptr or bool etc.
   * If cell doesn't exist default object is returned, which in case of a pointer is a nullptr.
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<UsesPtr>::type >
  const DataType* Get(const Vectorni& gidx) const {
    if (!Valid(gidx))
      return nullptr;

    return cells_[Id(gidx)].get();
  }


  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return Copy of contents of cell, could be shared_ptr or bool etc.
   * If cell doesn't exist default object is returned, which in case of a pointer is a nullptr.
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<!UsesPtr>::type >
  DataType Get(const Vectornd& x, bool checkValid = true) const {
    if (checkValid)
      if (!Valid(x))
        return DataType();

    return cells_[Id(x)];
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return const ref to contents of cell.
   * If cell holds pointer and is not allocate then a null reference is returned
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<!UsesPtr>::type >
  DataType Get(int id) const {
    if (id<0 || id >= nc_)
      return DataType();

    return cells_[id];
  }

  /**
   * Get the cell at a given grid index
   * @param gidx grid index
   * @return Copy of contents of cell, could be shared_ptr or bool etc.
   * If cell doesn't exist default object is returned, which in case of a pointer is a nullptr.
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<!UsesPtr>::type >
  DataType Get(const Vectorni& gidx) const {
    if (!Valid(gidx))
      return DataType();
    return cells_[Id(gidx)];
  }

  /**
   * Set the data corresponding to the position x.
   * @param x point
   * @param data
   * @return was able to set data or not
   */
  bool Set(const Vectornd& x, const DataType& data)  {
    int id = Id(x);
    if(id<0)
      return false;
    set_cells(id, data, cells_store_ptr_type);
    return true;
  }

  /**
   * Set the data corresponding to the grid index.
   * @param gidx grid index of a point on the grid
   * @param data
   * @return was able to set data or not
   */
  bool Set(const Vectorni& gidx, const DataType& data)  {
    int id = Id(gidx);
    if(id<0)
      return false;
    set_cells(id, data, cells_store_ptr_type);
    return true;
  }

  /**
   * Set the data at the given id
   * @param id a non-negative id
   * @param data the content of a cell
   * @return true if it was able to set the data
   */
  bool Set(int id, const DataType& data) {
    if (id<0 || id >= nc_)
      return false;
    set_cells(id, data, cells_store_ptr_type);
    return true;
  }

  /**
   * Provides an interface to loop over all cells fast while operating on id(array id) and gidx(grid index) of that cell.
   * for a grid of any dimension. It's is fast for grid of dimensionality 1,2,3,4,5.
   *
   * @param fun any function that operates on id and gidx
   */
  void LoopOver(std::function<void(int id, const Vectorni& gidx)> fun) const {
    int id=0;
    Vectorni gidx;
    switch(n_){
      case 1:
        for(int i0=0; i0<gs_[0]; i0++){
          gidx[0] = i0;
          fun(id,gidx);
          id++;
        }
        break;
      case 2:
        for(int i1=0; i1<gs_[1]; i1++){
          for(int i0=0; i0<gs_[0]; i0++){
            gidx[0] = i0; gidx[1] = i1;
            fun(id,gidx);
            id++;
          }
        }
        break;
      case 3:
        for(int i2=0; i2<gs_[2]; i2++){
          for(int i1=0; i1<gs_[1]; i1++){
            for(int i0=0; i0<gs_[0]; i0++){
              gidx[0] = i0; gidx[1] = i1; gidx[2] = i2;
              fun(id,gidx);
              id++;
            }
          }
        }
        break;
      case 4:
        for(int i3=0; i3<gs_[3]; i3++){
          for(int i2=0; i2<gs_[2]; i2++){
            for(int i1=0; i1<gs_[1]; i1++){
              for(int i0=0; i0<gs_[0]; i0++){
                gidx[0] = i0; gidx[1] = i1; gidx[2] = i2; gidx[3] = i3;
                fun(id,gidx);
                id++;
              }
            }
          }
        }
        break;
      case 5:
        for(int i4=0; i4<gs_[4]; i4++){
          for(int i3=0; i3<gs_[3]; i3++){
            for(int i2=0; i2<gs_[2]; i2++){
              for(int i1=0; i1<gs_[1]; i1++){
                for(int i0=0; i0<gs_[0]; i0++){
                  gidx[0] = i0; gidx[1] = i1; gidx[2] = i2; gidx[3] = i3;  gidx[4] = i4;
                  fun(id,gidx);
                  id++;
                }
              }
            }
          }
        }
        break;
      default:
        Vectornp1i gidxe = Vectornp1i::Zero();//e stands for extra element
        Vectornp1i  gse; gse << gs_,0; //0 to indicate iteration over all cells over
        int dim = 0;
        while (gidxe[n_]==0) {
          gidx = gidxe.head(n_);
          fun(id,gidx); //This function is called in loop
          id++; gidxe[0]++;
          while(gidxe[dim]==gse[dim]){
            gidxe[dim]=0;
            gidxe[++dim]++;
            if(gidxe[dim]!=gse[dim]) //here it is compared to last element of gse
              dim=0;
          }
        }
    }
  }

  /**
   * Converts a set of points in metric coordinates to grid coordinates
   * @param metric_coordinates metric coordinates
   * @param grid_coordinates For 2D grid, the Grid coordinates are (0,0) for cell with id=0
   */
  void ToGridCoordinates(const std::vector<Vectornd>& metric_coord, std::vector<Vectornd>* grid_coord ) const{
    grid_coord->resize(metric_coord.size());
    Vectornd point0;
    CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <metric_coord.size(); i++)
      (*grid_coord)[i] = (metric_coord[i] - point0).array()/cs_.array();
  }

  /**
   * Converts a set of points in metric coordinate to grid coordinates in place
   * @param coord In metric coordinates. After update they are in grid coordinates
   */
  void ToGridCoordinates(std::vector<Vectornd>* coord) const{

    Vectornd point0;
    Vectorni gidx0 = Vectorni::Zero();
    CellCenter(gidx0, &point0);
    for(size_t i=0; i <coord->size(); i++)
      (*coord)[i] = ((*coord)[i] - point0).array()/cs_.array();
  }

  /**
   * Converts a point in metric coordinate to grid coordinates
   * @param metric_coord metric coordinate
   * @param grid_coord grid coordinate
   */
  void ToGridCoordinates(const Vectornd& metric_coord, Vectornd* grid_coord) const{
    Vectornd point0;
    CellCenter(point0,Vectorni::Zero());
    *grid_coord = (metric_coord - point0).array()/cs_.array();
  }

  /**
   * Converts grid coordinates to metric coordinate
   * @param grid_coord grid coordinate
   * @param metric_coord metric coordinate
   */
  void FromGridCoordinates(const std::vector<Vectornd>& grid_coord, std::vector<Vectornd>* metric_coord) const{
    metric_coord->resize(grid_coord.size());
    Vectornd point0;
    CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <metric_coord->size(); i++)
      (*metric_coord)[i] = grid_coord[i].array()*cs_.array() + point0;
  }

  /**
   * Convert grid coordinate to metric coordinate in place
   * @param coord In grid coordinates and after update they are in metric
   */
  void FromGridCoordinates(std::vector<Vectornd>* coord) const{
    Vectornd point0;
    CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <coord->size(); i++)
      (*coord)[i] = (*coord)[i].array()*cs_.array() + point0;
  }

  /**
   * Convert a grid coordinate to metric coordinate
   * @param grid_coord
   * @param metric_coord
   */
  void FromGridCoordinates(const Vectornd& grid_coord, Vectornd* metric_coord) const{
    Vectornd point0;
    CellCenter(point0,Vectorni::Zero());
    *metric_coord = grid_coord.array()*cs_.array() + point0;
  }


  /**
   * Serialize data of the map class and save it in a binary file
   * @param filename
   */
  void Save(const std::string& filename) const{
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::ofstream fs(filename, std::fstream::out | std::ios::binary);
    if(fs.is_open()){

      //Copy the Grid data into the protobuff class
      dsl::ProtobufGrid pb;
      pb.set_n(n_);
      pb.set_nc(nc_);

      pb.mutable_xlb()->Reserve(n_);
      pb.mutable_xub()->Reserve(n_);
      pb.mutable_ds()->Reserve(n_);
      pb.mutable_cs()->Reserve(n_);
      pb.mutable_gs()->Reserve(n_);
      pb.mutable_cgs()->Reserve(n_);
      pb.mutable_wd()->Reserve(n_);
      for(int i=0; i < n_; i++){
        pb.mutable_xlb()->AddAlreadyReserved(xlb_[i]);
        pb.mutable_xub()->AddAlreadyReserved(xub_[i]);
        pb.mutable_ds()->AddAlreadyReserved(ds_[i]);
        pb.mutable_cs()->AddAlreadyReserved(cs_[i]);
        pb.mutable_gs()->AddAlreadyReserved(gs_[i]);
        pb.mutable_cgs()->AddAlreadyReserved(cgs_[i]);
        pb.mutable_wd()->AddAlreadyReserved(wd_[i]);
      }

      CellsToPb(pb, cells_store_ptr_type);

      pb.SerializeToOstream(&fs);

      fs.close();
      google::protobuf::ShutdownProtobufLibrary();
    }else{
      std::cout<<"Couldn't open file:"<< filename<<" to write"<<std::endl;
    }

  }

  /**
   * Read data from a binary file, deserialize data and create a map object from it
   * @param filename
   * @return The object loaded
   */
  static Ptr Load(const std::string& filename){
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    //google::protobuf::io::CodedInputStream::SetTotalBytesLimit(500000000,250000000);
    Ptr grid;
    std::ifstream fs (filename, std::fstream::in | std::ios::binary);
    if(fs.is_open()){
      dsl::ProtobufGrid pb;
      if(!pb.ParseFromIstream(&fs))
        return nullptr;
      fs.close();


      // size checks
      int n_allocated_cells;
      if(cells_store_ptr)
        n_allocated_cells = pb.ids_allocated_size();
      else
        n_allocated_cells = pb.nc();

      if( pb.data_size()  != n_allocated_cells){
        throw std::length_error("data field is of wrong size");
        return nullptr;
      }

      if(pb.xlb_size() != pb.n() || pb.xub_size()   != pb.n() ||
          pb.ds_size()  != pb.n() || pb.cs_size()    != pb.n() ||
          pb.gs_size()  != pb.n() || pb.cgs_size()   != pb.n() ||
          pb.wd_size()  != pb.n() ){
        throw std::length_error("All or one of xlb, xub, ds, cs, gs, cgs and wd is of wrong size");
        return nullptr;
      }

      // Data probably is valid, so initilize the grid and load data into it
      grid.reset(new Grid());

      if(grid->n_ !=  pb.n()){ //one last check
        throw std::length_error("pb's n doesn't match with grid.h");
        return nullptr;
      }

      grid->nc_ =  pb.nc();
      grid->cells_.resize(pb.nc());

      for(int i=0; i < grid->n_; i++){
        grid->xlb_[i] = pb.xlb(i);
        grid->xub_[i] = pb.xub(i);
        grid->ds_[i]  = pb.ds(i);
        grid->cs_[i]  = pb.cs(i);
        grid->gs_[i]  = pb.gs(i);
        grid->cgs_[i] = pb.cgs(i);
        grid->wd_[i]  = pb.wd(i);
      }
      PbToCells(*grid,pb,cells_store_ptr_type);


    }else{
      std::cout<<"Couldn't open file:"<< filename<<" to read"<<std::endl;
    }

    return grid;
  }


  // Accessors
  inline const int&       n(void)   const {
    return n_;
  }
  inline const int&       nc(void)  const {
    return nc_;
  }
  inline const Vectornd& xlb(void) const {
    return xlb_;
  }
  inline const Vectornd& xub(void) const {
    return xub_;
  }
  inline const Vectornd& ds(void) const {
    return ds_;
  }
  inline const Vectornd& cs(void)  const {
    return cs_;
  }
  inline const Vectorni& gs(void)  const {
    return gs_;
  }
  inline const Vectorni& cgs(void)  const {
    return cgs_;
  }
  inline const Vectornb& wd(void)  const {
    return wd_;
  }
  /**
   * Getter function for cells. This method is only enabled if cells don't store unique_ptr
   * @param id
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<!UsesPtr>::type >
  inline const std::vector<CellType>& cells(void)  const {
    return cells_;
  }

  // Mutators
  /**
   * setter function for cells
   * @param vals
   */
  inline void set_cells(const std::vector<DataType>& vals){
    if(nc_ != vals.size())
      throw std::length_error(" cells don't have the same length as nc.");
    set_cells(vals, cells_store_ptr_type);
  }

  /**
   * This method is only enabled if cells store shared_ptr
   * @param id
   */
  template< bool UsesPtr = UsePtr, typename = typename std::enable_if<UsesPtr>::type >
  void delete_cell(int id)
  {
    cells_.at(id).reset();
  }

  // To give data access to Slice and Stack
  template < class PointType_, class DataType_, bool UsePtr_> friend class Grid;

private:
  /**
   * Private constructor only to be used by the Load method
   */
  Grid(){}

  /**
   * setter function for cells_[id] if holds DataType and not pointer/smart_ptr to DataType
   * @param id id of the cell
   * @param val value
   * @param
   */
  inline void set_cells(int id, const DataType& val, std::false_type){
    cells_.at(id) = val;
  }

  /**
   * setter function for cells_[id] if holds pointer/smart_ptr to DataType
   * @param id id of the cell
   * @param val value
   * @param
   */
  inline void set_cells(int id, const DataType& val, std::true_type){
    cells_.at(id).reset(new DataType(val));
  }

  /**
   * setter function for cells if cells hold DataType and not pointer/smart_ptr to DataType
   * @param vals vector of values
   * @param
   */
  inline void set_cells(const std::vector<DataType>& vals, std::false_type){
    cells_ = vals;
  }

  /**
   * setter function for cells if it holds pointer/smart_ptr to DataType
   * @param vals vector of values
   * @param
   */
  inline void set_cells(const std::vector<DataType>& vals, std::true_type){
    for(int id = 0; id < nc_ ; id++){
      cells_.at(id).reset(new DataType(vals[id]));
    }
  }

  /**
   * Loads data from cells to protocol buffer for non shared_ptr CellType
   * @param pb
   * @param false_type
   */
  void CellsToPb(dsl::ProtobufGrid& pb, std::false_type) const {
    pb.mutable_data()->Reserve(nc_);
    for(int id=0; id < nc_; id++){
      std::string str;
      ValToString(cells_[id], &str, std::is_pod<DataType>());
      pb.add_data(str);
    }
  }

  /**
   * Loads data from cells to protocol buffer shared_ptr CellType
   * @param pb
   * @param true_type
   */
  void CellsToPb(dsl::ProtobufGrid& pb, std::true_type) const {
    if(!nc_)
      return;

    pb.mutable_ids_allocated()->Reserve(nc_); //reserve max
    pb.mutable_data()->Reserve(nc_); //reserve max
    for(int id=0; id < nc_; id++){
      if(cells_[id]){
        pb.add_ids_allocated(id);

        std::string str;
        ValToString(*cells_[id], &str, std::is_pod<DataType>());
        pb.add_data(str);
      }
    }
  }

  /**
   * converts a DataType to string for pod type
   * @param val
   * @param ss
   * @param
   */
  void ValToString(const DataType& val, std::string* str, std::true_type) const {
    int n_bytes = sizeof(DataType);
    str->resize(n_bytes);
    str->replace(0,n_bytes,(char*)&val, n_bytes);
  }

  /**
   * converts a DataType to string for non-pod type
   * @param val
   * @param ss
   * @param
   */
  void ValToString(const DataType& val, std::string* str, std::false_type) const{
    val.SerializeToString(str);
  }


  /**
   * Update cells_ from Protobuf for cell directly holding DataType
   * @param grid
   * @param pb
   * @param
   */
  static void PbToCells(Grid& grid, dsl::ProtobufGrid& pb, std::false_type){
    for(int id=0; id < grid.nc_; id++){
      DataType val;
      StringToVal(pb.data(id),&val, std::is_pod<DataType>());
      grid.cells_[id] = val;
    }
  }
  /**
   * Update cells_ from Protobuf for cell holding shared_ptr to DataType
   * @param grid
   * @param pb
   * @param
   */
  static void PbToCells(Grid& grid, dsl::ProtobufGrid& pb, std::true_type){
    for(int i=0; i < pb.ids_allocated_size(); i++){
      int id = pb.ids_allocated(i);
      DataType val;
      StringToVal(pb.data(i),&val, std::is_pod<DataType>());
      grid.cells_[id].reset(new DataType(val));
    }
  }


  /**
   * Converts a DataType data to string for pod type
   * @param ss
   * @param val
   * @param
   */
  static void StringToVal(const std::string& str, DataType* val, std::true_type){
    int n_bytes = sizeof(DataType);
    std::memcpy(val,str.c_str(), n_bytes);
  }

  /**
   * Converts a DataType data to string for non-pod type
   * @param ss
   * @param val
   * @param
   */
  static void StringToVal(const std::string& str, DataType* val, std::false_type){
    val->ParseFromString(str);
  }

  const int n_ = PointType::SizeAtCompileTime; ///< grid dimension
  int nc_ = 0;   ///< number of cells in grid
  Vectornd xlb_; ///< state lower bound
  Vectornd xub_; ///< state upper bound
  Vectornd ds_;  ///< dimension size (ds=xub-xlb)
  Vectornd cs_;  ///< cell length size per dimension
  Vectorni gs_;  ///< number of cells per dimension
  Vectorni cgs_; ///< cumulative(product) of gs. For n=3, cgs = [1, gs[0], gs[0]*gs[1]]
  Vectornb wd_;  ///< which dimensions are wrapped
  std::vector<CellType> cells_; ///< grid cells

  //const DataType default_val_ = DataType(); ///< default values
  const DataType& null_ref_ = *(DataType*)0; ///< null reference
};

/**
 * An n-dimenensional occupancy map storing type T in every cell. T=Bool
 * is works well to represent the occupancy
 */
template <typename T, int n>
using Map = Grid< Eigen::Matrix<double,n,1>, T, false >;

template < class PointType, class DataType>
using Lattice = Grid<PointType, DataType, true>;

}

#endif
