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
* traversibility of a cell
*
* Note that this data structure is only viable up to a few dimensions,
* e.g. dim=5 or 6.
* TODO: Make sure origin lies at center of a cell
 */
template < class PointType, class CellContent>
 class GridCore {
 public:

 using Vectornd =  Eigen::Matrix< double, PointType::SizeAtCompileTime, 1 >;     //   n  dim eigen vector of doubles
 using Vectorni =  Eigen::Matrix< int, PointType::SizeAtCompileTime, 1 >;        //   n  dim eigen vector of ints

 using Vectornm1d =  Eigen::Matrix< double, PointType::SizeAtCompileTime-1, 1 >; //(n-1) dim eigen vector of doubles
 using Vectornm1i =  Eigen::Matrix< int, PointType::SizeAtCompileTime-1, 1 >;    //(n-1) dim eigen vector of ints

 using Vectornp1d =  Eigen::Matrix< double, PointType::SizeAtCompileTime+1, 1 >; //(n+1) dim eigen vector of doubles
 using Vectornp1i =  Eigen::Matrix< int, PointType::SizeAtCompileTime+1, 1 >;    //(n+1) dim eigen vector of ints

 using GridCorePtr = shared_ptr< GridCore<PointType,CellContent> >;

 using GridCoreSlice = GridCore<Vectornm1d,CellContent>;
 using GridCoreSlicePtr = shared_ptr<GridCoreSlice>;

 using GridCoreExtradim = GridCore<Vectornp1d,CellContent>;
 using GridCoreExtradimPtr = shared_ptr<GridCoreExtradim>;

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
 Vectorni cgs; ///< cumulative(product) number of cells. For n=3, cgs = [1, gs[0], gs[0]*gs[1]]
 Vectorni wd;  ///< which dimensions are wrapped
 vector<CellContent> cells; //cells


  /**
   * Initialize the grid using state lower bound, state upper bound, the number of grid cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of grid cells per each dimension
   * @param wd indicates if a dimension if wrapped or not. wd[i]=0(flat dim) wd[i]>0(wrapped dim).
   * For wrapped dimension i, the ds[i](dimension size) is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   */
  GridCore(const Vectornd& xlb, const Vectornd& xub, const Vectorni& gs, const Vectorni wd = Vectorni::Zero())
      : n(xlb.size()), xlb(xlb), xub(xub), gs(gs), wd(wd), cells(0){
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
  }

  /**
   * Initialize the map using state lower bound, state upper bound, and a suggested cell size
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param scs suggested cell size along each dimension. The true cs is the one which results in
   * integer number of grid cells between xlb and xub
   * @param wd indicates if a dimension if wrapped or not. wd[i]=0(flat dim) wd[i]>0(wrapped dim).
   * For wrapped dimension i, the ds[i](dimension size) is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   */
  GridCore(const Vectornd& xlb, const Vectornd& xub, const Vectornd& scs, const Vectorni wd = Vectorni::Zero())
      : n(xlb.size()), xlb(xlb), xub(xub), wd(wd), cells(0){
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
  }

  /**
   * Initialize the map using a suggested state lower bound, suggested state upper bound, and a suggested cell size.
   *
   * @param sxlb state lower bound
   * @param sxub state upper bound
   * @param scs suggested cell size along each dimension. The true cs is the one which results in
   * integer number of grid cells between xlb and xub
   * @param wd indicates if a dimension if wrapped or not. wd[i]=0(flat dim) wd[i]>0(wrapped dim).
   * For wrapped dimension i, the ds[i](dimension size) is given by xub[i]-xlb[i] (e.g. for angles it ds[i]=2*M_PI)
   */
  GridCore(const Vectornd& sxlb, const Vectornd& sxub, const Vectornd& scs, const Vectorni wd,const GridConfig& cfg)
      : n(sxlb.size()), wd(wd), cells(0){
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
  }

  GridCore(const GridCore &gridcore)
  : n(gridcore.n), nc(gridcore.nc), xlb(gridcore.xlb), xub(gridcore.xub), ds(gridcore.ds),
    cs(gridcore.cs), gs(gridcore.gs), cgs(gridcore.cgs), wd(gridcore.wd) {
    cells.resize(nc);
    //don't want to copy smart pointers
    //if(std::is_copy_constructible<CellContent>::value) //TODO: Should this be used?
    if(std::is_arithmetic<CellContent>::value)
      cells = gridcore.cells;
  }

  virtual ~GridCore() {}

//  /**
//   * Adds an extra dimension to the grid and repeat all other dimensions
//   * @param dim dim coordinate/dimension index
//   * @param lb lower bound along that dimension
//   * @param ub upper bound along that dimension
//   * @param wd dimension is wrapped or not
//   * @return
//   */
//  GridCoreExtradimPtr addDim(int dim, double lb, double ub, bool wrapped){
//
//  }


  /**
   * Create a slice of the current grid along dimension at particular index
   * @param idx index in a grid along the dim coordinate/dimension
   * @param dim coordinate/dimension index
   * @return
   */
  GridCoreSlicePtr Slice(int idx, int dim) const {
    Vectornm1d xlb_slice, xub_slice;
    Vectornm1i gs_slice,wd_slice;

    //Get the new dimension
    xlb_slice = removeDim(xlb,dim);
    xub_slice = removeDim(xub,dim);
    gs_slice  = removeDim(gs,dim);
    wd_slice  = removeDim(wd,dim);

    GridCoreSlicePtr pslice(new GridCoreSlice(xlb_slice,xub_slice, gs_slice, wd_slice));

    for(int id_slice=0; id_slice < pslice->nc; id_slice++){
      Vectornm1i midx_slice; pslice->Index(midx_slice,id_slice);
      Vectorni midx = insertDim(midx_slice,dim,idx); //midx is multidim index
      pslice->cells[id_slice] = Get(Id(midx));
    }
    return pslice;
  }

  /**
    * Create a slice of the current grid along dimension at particular index
    * @param idx index in a grid along the dim coordinate/dimension
    * @param dim coordinate/dimension index
    * @return
    */
   GridCoreSlicePtr Slice(double val, int dim) const {
     int idx = Index(val,dim);
     return Slice(idx,dim);
   }

  /**
   * Create new GridCore object whose cell resolution is scale times higher.
   * Useful for increasing resolution of a map for display
   * @param scale
   * @return
   */
  GridCorePtr ScaleUp(int scale) const {
    assert(scale>1);
    GridCorePtr pscaled;
    if(scale<=1)
      return pscaled;
    Vectorni gs_scaled = gs*scale;
    pscaled.reset(new GridCore<PointType, CellContent>(xlb,xub, gs_scaled, wd));

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
 * @param idx multidimensional index
 * @return true if within bounds
 */
  bool Valid(const Vectorni& idx) const{
    for(size_t i=0; i< n; i++){
      if(idx[i]<0 || idx[i]>=gs[i])
        return false;
    }
    return true;
  }

/**
 * Finds the CellCenter corresponding to the input idx
 * @param x The center point
 * @param idx Multidimensional index
 * @param checkValid Check if the index is within bounds
 * @return
 */
  bool CellCenter(Vectornd& x,const Vectorni& idx, bool checkValid = true) const{
    if(checkValid)
      if(!Valid(idx))
        return false;

    x = xlb.array() +  (idx.template cast<double>() + Vectornd::Constant(0.5)).array()*cs.array() ;
    return true;
  }

  /**
   * Finds the CellCenter corresponding to the input idx
   * @param x The updated center
   * @param id id of point in grid array
   * @param checkValid Check if the index is within bounds
   * @return
   */
  bool CellCenter(Vectornd& x, int id) const{
    if(id>=nc || id<0)
      return false;
    Vectorni idx;
    if(Index(idx,id)){
      x = xlb.array() +  (idx.template cast<double>()+Vectornd::Constant(0.5)).array()*cs.array() ;
      return true;
    }else{
      return false;
    }
  }

  /**
   * Returns cc[dim] where cc is the cell center of any cell with grid index = idx along dimension dim
   * @param idx grid index along dim i. To clarify 0<= idx < gs[i]
   * @param dim the Dimension number
   * @return Returns cc[dim]
   */
  double CellCenterIth(int idx, int dim) const{
    return xlb[dim] + (idx +0.5)*cs[dim];
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
   * @param idx index in a grid along all coordinate/dimension
   * @return a computed id
   */
  int Id(const Vectorni& idx) const {
    return idx.transpose()*cgs;
  }

  /**
   * Get the grid index of the point x along i coordinate index
   * @param x point
   * @param i coordinate/dimension index
   * @return index in a grid along i coordinate/dimension
   */
  int Index(const Vectornd& x, int i) const {
    if(!wd[i]){ //dimension is flat
      return floor((x[i] - xlb[i]) / ds[i] * gs[i]);
    }else{  //dimension is wrapped
      double xi = x[i];
      while(xi < xlb[i]){xi += ds[i];}
      while(xi > xub[i]){xi -= ds[i];}
      return floor((xi - xlb[i]) / ds[i] * gs[i]);
    }
  }

  /**
   * Get the grid index of i coordinate/dimension
   * @param xi value of a point at i coordinate/dimension
   * @param i coordinate/dimension index
   * @return index in a grid along i coordinate/dimension
   */
  int Index(double xi, int i) const {
    if(!wd[i]){ //dimension is flat
      return floor((xi - xlb[i]) / ds[i] * gs[i]);
    }else{  //dimension is wrapped
      while(xi < xlb[i]){xi += ds[i];}
      while(xi > xub[i]){xi -= ds[i];}
      return floor((xi - xlb[i]) / ds[i] * gs[i]);
    }
  }

  /**
   * Get the grid index along all the dimensions
   * @param idx index in a grid along all coordinate/dimension
   * @param x point
   */
  void Index(Vectorni& idx, const Vectornd& x) const {
    for(size_t i=0;i<n;i++){
      if(!wd[i]){ //dimension is flat
        idx[i] = floor((x[i] - xlb[i]) / ds[i] * gs[i]);
      }else{  //dimension is wrapped
        double xi = x[i];
        while(xi < xlb[i]){xi += ds[i];}
        while(xi > xub[i]){xi -= ds[i];}
        idx[i] = floor((xi - xlb[i]) / ds[i] * gs[i]);
      }
    }
  }

  /**
   * Get the grid index along all the dimensions from the direct lookup id
   * @param id id for direct lookup in the grid array
   * @return false if index is out of range
   */
  bool Index(Vectorni& idx, int id) const {
    if(id>=nc || id<0)
      return false;
    for(int i=n-1; i>=0; i--){
      idx[i] = floor(id/cgs[i]);
      id -= idx[i]*cgs[i];
    }
    return true;
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
   static GridCorePtr Load(const string& filename) {
     GridCorePtr pgrid;
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

   friend class cereal::access;//give cereal access to non-public methods

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
