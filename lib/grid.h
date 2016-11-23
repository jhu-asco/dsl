// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu> and Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID_H
#define DSL_GRID_H

#include "cell.h"
#include "utils.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <fstream>

#include "cereal/archives/binary.hpp" //serialization support
#include "cereal/types/vector.hpp" //type support
#include "cereal/types/eigen.hpp" //type support

namespace dsl {

struct EmptyData {};
using std::shared_ptr;
using std::vector;
using namespace Eigen;

/**
 * An n-dimenensional grid consisting of abstract "cells", or elements
 * identified by a set of coordinates of type PointType, each cell
 * containing data of type CellContent. The CellContent can directly be
 * the intended data or shared_ptr to intended data. Using the shared_ptr
 * is useful to avoid memory allocation until it's actually needed.
 * CellContent=bool/double can be used to make a grid of occupancy or
 * traversibility of a cell.
 *
 * When fastindex is set to true then GridCore class stores the grid
 * location for every cell. This is to provide fast conversion of array id
 * to grid index. For 1000x1000x100 grid, Index() is 2 orders of magnitude
 * faster at the cost of using additional memory
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 *
 * TODO: 1. Make sure origin lies at center of a cell
 *       2. Check for bound in GetSlice and GetStack methods
 *       3. Implement GetStack method that can perform occupancy dilation
 *
 */
template < class PointType, class CellContent>
class GridCore {
public:


  //n-dimensional vectors
  using Vectornd   =  Eigen::Matrix< double,  PointType::SizeAtCompileTime,   1 >;
  using Vectorni   =  Eigen::Matrix< int,     PointType::SizeAtCompileTime,   1 >;
  using Vectornu   =  Eigen::Matrix< uint16_t,PointType::SizeAtCompileTime,   1 >;
  using Vectornb   =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime,   1 >;

  //n-1 dimensional vectors
  using Vectornm1d =  Eigen::Matrix< double,  PointType::SizeAtCompileTime-1, 1 >;
  using Vectornm1i =  Eigen::Matrix< int,     PointType::SizeAtCompileTime-1, 1 >;
  using Vectornm1u =  Eigen::Matrix< uint16_t,PointType::SizeAtCompileTime-1, 1 >;
  using Vectornm1b =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime-1, 1 >;

  //n+1 dimensional vectors
  using Vectornp1d =  Eigen::Matrix< double,  PointType::SizeAtCompileTime+1, 1 >;
  using Vectornp1i =  Eigen::Matrix< int,     PointType::SizeAtCompileTime+1, 1 >;
  using Vectornp1u =  Eigen::Matrix< uint16_t,PointType::SizeAtCompileTime+1, 1 >;
  using Vectornp1b =  Eigen::Matrix< bool,    PointType::SizeAtCompileTime+1, 1 >;

  using Ptr = shared_ptr< GridCore<PointType,CellContent> >;

  using Slice = GridCore<Vectornm1d,CellContent>;
  using SlicePtr = shared_ptr<Slice>;

  using Stack = GridCore<Vectornp1d,CellContent>;
  using StackPtr = shared_ptr<Stack>;

  /**
   * Configuration to control the grid sizes, cell sizes and grid bounds.
   */
  enum GridConfig{
    ORG_AT_CELL_CENTER, ///< makes sure the origin lies in the center of a gridcell
  };


  int n;        ///< grid dimension
  int nc = 0;   ///< number of cells in grid
  Vectornd xlb; ///< state lower bound
  Vectornd xub; ///< state upper bound
  Vectornd ds;  ///< dimension size (ds=xub-xlb)
  Vectornd cs;  ///< cell length size per dimension
  Vectorni gs;  ///< number of cells per dimension
  Vectorni cgs; ///< cumulative(product) of gs. For n=3, cgs = [1, gs[0], gs[0]*gs[1]]
  Vectornb wd;  ///< which dimensions are wrapped
  bool fi; ///< enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
  vector<Vectornu> indices; ///< grid index of all the cells
  vector<CellContent> cells; ///< grid cells

  /**
   * Initialize the grid using state lower bound, state upper bound, the number of grid cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of grid cells per each dimension
   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   */
  GridCore(const Vectornd& xlb, const Vectornd& xub, const Vectorni& gs, const Vectornb wd = Vectornb::Zero(), bool fi = false)
  : n(xlb.size()), xlb(xlb), xub(xub), gs(gs), wd(wd),fi(fi),indices(0), cells(0){
    ds = xub - xlb;
    nc = 1;
    cgs[0] = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      nc *= gs[i]; // total number of cells
      cs[i] = ds[i] / gs[i];
      if(i>0)
        cgs[i] = cgs[i-1]*gs[i-1];
    }
    cells.resize(nc);
    if(fi)
      ResetIndices();
  }

  /**
   * Initialize the map using state lower bound, state upper bound, and a suggested cell size
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param scs suggested cell size along each dimension. The true cs is the one which results in
   * integer number of grid cells between xlb and xub
   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   */
  GridCore(const Vectornd& xlb, const Vectornd& xub, const Vectornd& scs, const Vectornb wd = Vectornb::Zero(),bool fi = false)
  : n(xlb.size()), xlb(xlb), xub(xub), wd(wd),fi(fi),indices(0), cells(0){
    ds = xub - xlb;
    nc = 1;
    cgs[0] = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] < xub[i]);
      assert(scs[i] > 0);
      gs[i] = floor(ds[i] / scs[i]);
      cs[i] = ds[i] / gs[i];
      nc *= gs[i]; // total number of cells
      if(i>0)
        cgs[i] = cgs[i-1]*gs[i-1];
    }
    cells.resize(nc);
    if(fi)
      ResetIndices();
  }

  /**
   * Initialize the map using a suggested state lower bound, suggested state upper bound, and a suggested cell size.
   *
   * @param sxlb state lower bound
   * @param sxub state upper bound
   * @param scs suggested cell size along each dimension. The true cs is the one which results in
   * integer number of grid cells between xlb and xub
   * @param wd indicates if a dimension if wrapped or not. wd[i]=false(flat dim) wd[i]=1(wrapped dim).
   * For wrapped dimension i, the ds[i],i.e. dimension size, is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   */
  GridCore(const Vectornd& sxlb, const Vectornd& sxub, const Vectornd& scs, const Vectornb wd,const GridConfig& cfg,bool fi = false)
  : n(sxlb.size()), wd(wd),fi(fi),indices(0), cells(0){
    nc = 1;
    cgs[0] = 1;
    for (int i = 0; i < n; ++i) {
      assert(sxlb[i] < sxub[i]);
      assert(scs[i] > 0);

      switch(cfg){
        case ORG_AT_CELL_CENTER:
          if(!wd[i]){ //if dimension is flat(cs=scs)
            cs[i] = scs[i];
            int n_org2ub = floor(sxub[i]/cs[i] - 0.5);
            int n_org2lb =  ceil(sxlb[i]/cs[i] + 0.5);
            gs[i] = n_org2ub  - n_org2lb + 1;
            xub[i] = (n_org2ub +0.5)*cs[i];
            xlb[i] = (n_org2lb -0.5)*cs[i];
            ds[i] = xub[i] - xlb[i];
          }else{ //dimension is wrapped
            ds[i] = sxub[i] - sxlb[i]; //For angles this is = 2*M_PI
            int gs[i] = 2*ceil(ds[i]/(2*scs[i])); //Even number of cells
            cs[i] = ds[i]/gs[i];
            xlb[i] = -ds[i]/2 + cs[i]/2; //shifted by half cell so that pi/(pi/2)/0 lies at center of a cell
            xub[i] =  ds[i]/2 + cs[i]/2; //shifted by half cell so that pi/(pi/2)/0 lies at center of a cell
          }
          break;
        default: // Defaults to what's done in other constructors
          assert("Not implemented yet" && 0);
          xub[i] = sxub[i];
          xub[i] = sxub[i];
          ds[i] = xub[i] - xlb[i];
          gs[i] = floor(ds[i] / scs[i]);
          cs[i] = ds[i] / gs[i];
          break;
      }

      nc *= gs[i]; // total number of cells
      if(i>0)
        cgs[i] = cgs[i-1]*gs[i-1];
    }
    cells.resize(nc);
    if(fi)
      ResetIndices();
  }

  /**
   * Copies the grid structure and allocates memory for all but assigns only for arithmetic CellContent
   * @param gridcore
   */
  GridCore(const GridCore &gridcore)
  : n(gridcore.n), nc(gridcore.nc), xlb(gridcore.xlb), xub(gridcore.xub), ds(gridcore.ds),
    cs(gridcore.cs), gs(gridcore.gs), cgs(gridcore.cgs), wd(gridcore.wd),
    fi(gridcore.fi), indices(gridcore.indices){
    cells.resize(nc);

    //don't copy smart pointers
    if(std::is_arithmetic<CellContent>::value) //TODO: should is_copy_constructible<CellContent>::value be used?
      cells = gridcore.cells;
  }

  virtual ~GridCore() {}

  /**
   * Copies the grid structure and allocates memory for all but assigns only for arithmetic CellContent
   * @param The grid to copy the structure from
   * @return
   */
  GridCore& operator=(const GridCore& gridcore){
    n   = gridcore.n;
    nc  = gridcore.nc;
    xlb = gridcore.xlb;
    xub = gridcore.xub;
    gs  = gridcore.gs;
    cs  = gridcore.cs;
    cgs = gridcore.cgs;
    wd  = gridcore.wd;
    fi = gridcore.fi;
    indices = gridcore.indices;
    cells.resize(nc);
    if(std::is_arithmetic<CellContent>::value) //TODO: should is_copy_constructible<CellContent>::value be used?
      cells = gridcore.cells;
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
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   * @param init initialize the Stack from grid
   * @return
   */
  StackPtr GetStack(int i, double lbi, double ubi, int gsi, bool wdi = false, bool fi = false, bool init = true){
    Vectornp1d xlb_stack, xub_stack;
    Vectornp1i gs_stack;
    Vectornp1b wd_stack;

    //Get the new dimension
    xlb_stack = insertDim(xlb,i, lbi);
    xub_stack = insertDim(xub,i, ubi);
    gs_stack  = insertDim(gs, i, gsi);
    wd_stack  = insertDim(wd, i, wdi);

    StackPtr pstack(new Stack(xlb_stack,xub_stack, gs_stack, wd_stack,fi));

    if(!init)
      return pstack;

    for(int id=0; id < nc; id++){
      Vectorni gidx; Index(gidx,id); //grid index
      CellContent val = cells[id];
      for( int idx=0; idx < pstack->gs[i] ; idx++ ){
        Vectornp1i gidx_stack = insertDim(gidx,i,idx);
        pstack->cells[pstack->Id(gidx_stack)] = val;
      }
    }
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
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   * @param init initialize the Stack from grid
   * @return
   */
  StackPtr GetStack(int i, double lbi, double ubi, double scsi, bool wdi = false, bool fi = false, bool init = true){
    Vectornp1d xlb_stack, xub_stack,scs_stack;
    Vectornp1b wd_stack;


    //Get the new dimension
    xlb_stack = insertDim(xlb,i, lbi);
    xub_stack = insertDim(xub,i, ubi);
    scs_stack = insertDim(cs, i, scsi);
    wd_stack  = insertDim(wd, i, wdi);

    StackPtr pstack(new Stack(xlb_stack,xub_stack, scs_stack, wd_stack,fi));

    if(!init)
      return pstack;

    for(int id=0; id < nc; id++){
      Vectorni gidx; Index(gidx,id);
      CellContent val = cells[id];
      for( int idx=0; idx < pstack->gs[i] ; idx++ ){
        Vectornp1i gidx_stack = insertDim(gidx,i,idx);
        pstack->cells[pstack->Id(gidx_stack)] = val;
      }
    }
    return pstack;
  }

  /**
   * Create a slice of the current grid along dimension at particular index
   * @param idx The dimension along which the slice is created
   * @param i coordinated index/ dimension
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   * @param init initialize the slice from grid
   * @return a shared_ptr to the grid slice created
   */
  SlicePtr GetSlice(int idx, int i, bool fi = false, bool init = true) const {
    Vectornm1d xlb_slice, xub_slice;
    Vectornm1i gs_slice;
    Vectornm1b wd_slice;

    //Get the new dimension
    xlb_slice = removeDim(xlb,i);
    xub_slice = removeDim(xub,i);
    gs_slice  = removeDim(gs,i);
    wd_slice  = removeDim(wd,i);

    SlicePtr pslice(new Slice(xlb_slice,xub_slice, gs_slice, wd_slice,fi));

    if(!init)
      return pslice;

    for(int id_slice=0; id_slice < pslice->nc; id_slice++){
      Vectornm1i gidx_slice; pslice->Index(gidx_slice,id_slice);
      Vectorni gidx = insertDim(gidx_slice,i,idx);
      pslice->cells[id_slice] = Get(Id(gidx));
    }
    return pslice;
  }

  /**
   * Create a slice of the current grid along dimension at particular index
   * @param val val = point[i] where i is the dimension along which the slice is created
   * @param i coordinated index/ dimension
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   * @param init initialize the slice from grid
   * @return a shared_ptr to the grid slice created
   */
  SlicePtr GetSlice(double val, int i, bool fi = false, bool init = true) const {
    int idx = Index(val,i);
    return GetSlice(idx,i,fi,init);
  }

    /**
      * Create a slice of the current grid along dimension at particular index
      * @param slice reference to the slice
      * @param idx index in a grid along the dim coordinate/dimension
      * @param dim coordinated index/ dimension
      */
     void GetSlice(Slice& slice, int idx, int dim, bool init = true) const {
       //Get the new dimension
       slice.xlb = removeDim(xlb,dim);
       slice.xub = removeDim(xub,dim);
       slice.gs  = removeDim(gs,dim);
       slice.wd  = removeDim(wd,dim);
       slice.cs  = removeDim(cs,dim);
       slice.ds  = removeDim(ds,dim);
       slice.nc=1;
       slice.cgs[0]=1;
       for(int i = 0; i < n-1; i++){
         slice.nc *=slice.gs[i];
         if(i>0)
           slice.cgs[i] = slice.cgs[i-1]*slice.gs[i-1];
       }
       slice.cells.resize(slice.nc);//if data is already allocated resize does nothing

       if(slice.fi)
         slice.ResetIndices();

       if(!init)
         return;

       for(int id_slice=0; id_slice < slice.nc; id_slice++){
         Vectornm1i midx_slice; slice.Index(midx_slice,id_slice);
         Vectorni midx = insertDim(midx_slice,dim,idx); //midx is multidim index
         slice.cells[id_slice] = Get(Id(midx));
       }
     }


  /**
   * Create new GridCore object whose cell resolution is scale times higher.
   * Useful for increasing resolution of a map for display
   * @param scale
   * @param fi enable fast conversion of array id to grid index(2 orders of magnitude times faster but takes memory)
   * @param init initialize the scaled map from grid
   * @return
   */
  Ptr ScaleUp(int scale, bool fi = false, bool init = true) const {
    assert(scale>1);
    Ptr pscaled;
    if(scale<=1)
      return pscaled;
    Vectorni gs_scaled = gs*scale;
    pscaled.reset(new GridCore<PointType, CellContent>(xlb,xub, gs_scaled, wd,fi));

    if(!init)
      return pscaled;

    for(int id_scaled = 0; id_scaled <pscaled->nc; id_scaled++){
      Vectornd cc; pscaled->CellCenter(cc,id_scaled);
      pscaled->cells[id_scaled] = Get(cc,false);
    }
    return pscaled;
  }

  /**
   * Check if point x is within grid bounds
   * @param x point
   * @return true if within bounds
   */
  bool Valid(const Vectornd& x, double eps = 1e-10) const {
    for (int i = 0; i < x.size(); ++i) {
      if(!wd[i]){ //if dimension is flat
        if (x[i] < xlb[i] + eps)
          return false;
        if (x[i] > xub[i] - eps)
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
    for(size_t i=0; i< n; i++){
      if(gidx[i]<0 || gidx[i]>=gs[i])
        return false;
    }
    return true;
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param x The center point
   * @param gidx Multidimensional index
   * @param checkValid Check if the index is within bounds
   * @return
   */
  bool CellCenter(Vectornd& x,const Vectorni& gidx, bool checkValid = true) const{
    if(checkValid)
      if(!Valid(gidx))
        return false;

    x = xlb.array() +  (gidx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs.array() ;
    return true;
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param x The updated center
   * @param id id of point in the storage array
   * @param checkValid Check if the index is within bounds
   * @return
   */
  bool CellCenter(Vectornd& x, int id) const{
    if(id>=nc || id<0)
      return false;
    Vectorni gidx;
    if(Index(gidx,id)){
      x = xlb.array() +  (gidx.template cast<double>()+Vectornd::Constant(0.5)).array()*cs.array() ;
      return true;
    }else{
      return false;
    }
  }

  /**
   * Returns cc[i] where cc is the cell center of any cell with grid index = idx along dimension dim
   * @param idx grid index along dimension i. To clarify 0<= idx < gs[i]
   * @param i coordinated index/ dimension
   * @return Returns cc[dim]
   */
  double CellCenterIth(int idx, int i) const{
    return xlb[i] + (idx +0.5)*cs[i];
  }

  /**
   * Get an id of point x useful for direct lookup in the grid array
   * @param x point
   * @return a computed id
   */
  int Id(const Vectornd& x) const {
    Vectorni idx; Index(idx,x);
    return idx.transpose()*cgs;
  }

  /**
   * Get the id of a point corresponding to the multidimensional index
   * @param gidx grid index
   * @return a computed id
   */
  int Id(const Vectorni& gidx) const {
    return gidx.transpose()*cgs;
  }

  /**
   * Get the grid index of the point x along i coordinate index
   * @param x point
   * @param i coordinated index/ dimension
   * @return grid index along i coordinate/dimension
   */
  int Index(const Vectornd& x, int i) const {
    double xi = x[i];
    if(wd[i]){ //dimension is wrapped
      while(xi < xlb[i]){xi += ds[i];}
      while(xi > xub[i]){xi -= ds[i];}
    }
    return floor((xi - xlb[i]) / ds[i] * gs[i]);
  }

  /**
   * Get the grid index of i coordinate/dimension
   * @param xi value of a point at i coordinate/dimension
   * @param i coordinated index/ dimension
   * @return grid index along i coordinate/dimension
   */
  int Index(double xi, int i) const {
    if(wd[i]){ //dimension is wrapped
      while(xi < xlb[i]){xi += ds[i];}
      while(xi > xub[i]){xi -= ds[i];}
    }
    return floor((xi - xlb[i]) / ds[i] * gs[i]);
  }

  /**
   * Get the grid index along all the dimensions
   * @param gidx grid index
   * @param x point
   */
  void Index(Vectorni& gidx, const Vectornd& x) const {
    for(size_t i=0;i<n;i++){
      double xi = x[i];
      if(wd[i]){ //dimension is wrapped
        while(xi < xlb[i]){xi += ds[i];}
        while(xi > xub[i]){xi -= ds[i];}
      }
      gidx[i] = floor((xi - xlb[i]) / ds[i] * gs[i]);
    }
  }

  /**
   * Get the grid index along all the dimensions from the direct lookup id
   * @param gidx grid index
   * @param id id for direct lookup in the grid array
   * @return false if index is out of range
   */
  bool Index(Vectorni& gidx, int id) const {
    if(id>=nc || id<0)
      return false;
    if(fi){
      gidx = indices[id].template cast<int>();
    }else{
      for(int i=n-1; i>=0; i--){
        gidx[i] = int(id/cgs[i]);
        id -= gidx[i]*cgs[i];
      }
    }
    return true;
  }

  /**
   * Saves the grid index for every cell to provide fast conversion from
   * array id to grid index
   */
  void ResetIndices(){
    indices.resize(nc);
    //Iterating over all cells starting from zero index
    //e stands for extra element
    Vectornu gsu = gs.template cast<uint16_t>();
    Vectornp1u idxe = Vectornp1u::Zero();
    Vectornp1u  gse; gse << gsu,1; //1 to indicate iteration over
    int i = 0;
    int id = 0;
    while (idxe[n]==0) {
      indices[id++] = idxe.head(n);
      idxe[0]++;
      while(idxe[i]==gse[i]){
        idxe[i]=0;
        idxe[++i]++;
        if(idxe[i]!=gse[i])
          i=0;
      }
    }
  }

  /**
   * Converts a set of points in metric coordinates to grid coordinates
   * @param grid_coordinates For 2D grid, the Grid coordinates are (0,0) for cell with id=0
   * @param metric_coordinates metric coordinates
   */
  void ToGridCoordinates(vector<Vectornd>& grid_coord, const vector<Vectornd>& metric_coord) const{
    grid_coord.resize(metric_coord.size());
    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <metric_coord.size(); i++)
      grid_coord[i] = (metric_coord[i] - point0).array()/cs.array();
  }

  /**
   * Converts a set of points in metric coordinate to grid coordinates in placeS
   * @param coord In metric coordinates. After update they are in grid coordinates
   */
  void ToGridCoordinates(vector<Vectornd>& coord) const{

    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <coord.size(); i++)
      coord[i] = (coord[i] - point0).array()/cs.array();
  }

  /**
   * Converts a point in metric coordinate to grid coordinates
   * @param grid_coord grid coordinate
   * @param metric_coord metric coordinate
   */
  void ToGridCoordinates(Vectornd& grid_coord, const Vectornd& metric_coord) const{
    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    grid_coord = (metric_coord - point0).array()/cs.array();
  }

  /**
   * Converts grid coordinates to metric coordinate
   * @param metric_coord
   * @param grid_coord
   */
  void FromGridCoordinates(vector<Vectornd>& metric_coord, const vector<Vectornd>& grid_coord) const{
    metric_coord.resize(grid_coord.size());
    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <metric_coord.size(); i++)
      metric_coord[i] = grid_coord[i].array()*cs.array() + point0;
  }

  /**
   * Convert grid coordinate to metric coordinate in place
   * @param coord In grid coordinates and after update they are in metric
   */
  void FromGridCoordinates(vector<Vectornd>& coord) const{
    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    for(size_t i=0; i <coord.size(); i++)
      coord[i] = coord[i].array()*cs.array() + point0;
  }

  /**
   * Convert a grid coordinate to metric coordinate
   * @param metric_coord
   * @param grid_coord
   */
  void FromGridCoordinates(Vectornd& metric_coord, const Vectornd& grid_coord) const{
    Vectornd point0; CellCenter(point0,Vectorni::Zero());
    metric_coord = grid_coord.array()*cs.array() + point0;
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return pointer to a cell or 0 if cell does not exist
   */
  CellContent Get(const Vectornd& x, bool checkValid = true) const {
    if (checkValid)
      if (!Valid(x))
        return 0;
    return Get(Id(x));
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  CellContent Get(int id) const {
    assert(id >= 0);
    if (id >= nc)
      return 0;
    return cells[id];
  }

  /**
   * Set the data corresponding to the position x.
   * @param x point
   * @param data
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   */
  bool Set(const Vectornd& x, const CellContent& data, bool checkValid = true)  {
    if (checkValid)
      if (!Valid(x))
        return false;

    int id = Id(x);
    assert(id >= 0 && id < nc);
    cells[id] = data;
    return true;
  }

  /**
   * Set the data corresponding to the multidimensional index.
   * @param idx multidimensional index of a point on the grid
   * @param data
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   */
  bool Set(const Vectorni& idx, const CellContent& data, bool checkValid = true)  {
    if (checkValid)
      if (!Valid(idx))
        return false;

    int id = IndexToId(idx);
    assert(id >= 0 && id < nc);
    cells[id] = data;
    return true;
  }

  /**
   * Set the data at the given id
   * @param id a non-negative id
   * @param data the content of a cell
   * @return true if it was able to set the data
   */
  bool Set(int id, const CellContent& data) {
    assert(id >= 0 && id<nc);
    if (id<0 || id >= nc)
      return false;
    cells[id] = data;
    return true;
  }

  /**
   * Serialize data of the map class and save it in a binary file
   * @param map
   * @param filename
   */
  static void Save( GridCore<Vectornd,CellContent>& grid, const string& filename) {
    std::ofstream fs(filename, std::fstream::out | std::ios::binary);
    assert(fs.is_open());
    if(fs.is_open()){
      //Serialize and write to file
      cereal::BinaryOutputArchive cerealout(fs);
      cerealout(grid);
      fs.close();
    }else{
      cout<<"Couldn't open file:"<< filename<<" to write"<<endl;
    }
  }


  /**
   * Read data from a binary file, deserialize data and create a map object from it
   * @param filename
   * @return
   */
  static Ptr Load(const string& filename) {
    Ptr pgrid;
    std::ifstream fs (filename, std::fstream::in | std::ios::binary);
    assert(fs.is_open());
    if(fs.is_open()){
      pgrid.reset(new GridCore<Vectornd,CellContent>());
      cereal::BinaryInputArchive cerealin(fs);
      cerealin(*pgrid);
    }else{
      cout<<"Couldn't open file:"<< filename<<" to read"<<endl;
    }
    return pgrid;
  }

private:
  /**
   * Private constructor only to be used by the Cereal serialization
   */
  GridCore(){}

  /**
   * give cereal access to non-public methods
   *
   */
  friend class cereal::access;

  /**
   * serealizes data members for cereal library
   * @param ar cereal(library) archive
   */
  template <class Archive>
  void serialize( Archive & ar ){
    ar( n, nc, xlb, xub, ds, cs, gs, cgs, wd, cells);
  }

};

/**
 * An n-dimenensional grid consisting of abstract "cells", or elements
 * identified by a set of coordinates of type PointType, each cell
 * containing data of type DataType.
 * A grid provides instant access to the elements by maintaining
 * an n-dimensional array of pointers to cells. Cells that are empty,
 * e.g. that are inside obstacles, correspond to null pointers.
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 */
template < class PointType, class DataType = EmptyData>
using Grid = GridCore< PointType,shared_ptr< Cell<PointType, DataType> > >;


/**
 * An n-dimenensional occupancy map storing type T in every cell. T=Bool
 * is works well to represent the occupancy
 */
template <typename T, int n>
using Map = GridCore< Matrix<double,n,1>, T >;



}

#endif
