// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRID_H
#define DSL_GRID_H

#include "cell.h"

namespace dsl {

  using namespace Eigen;  
  
  template<int n>
    class Grid {
    
    typedef Matrix<double, n, 1> Vectornd;
    typedef Matrix<int, n, 1> Vectorni;
    
  public:
    
  Grid(const Vectornd &xlb, const Vectornd &xub, const Vectorni &gs) :
    xlb(xlb), xub(xub), gs(gs) {
      nc = 1;
      for (int i = 0; i < n; ++i) {
        assert(xlb[i] <= xub[i]);
        assert(gs[i] > 0);
        nc *= gs[i];  // total number of cells
        cs[i] = (xub[i] - xlb[i])/gs[i];
      }

      //      cs = (xub - xlb).cwiseQuotient(gs);

      cells = new Cell<n>*[nc];
      memset(cells, 0, nc*sizeof(Cell<n>*)); // initialize all of them nil
    }

    virtual ~Grid()  {
      for (int i = 0; i < nc; ++i)
        delete cells[i];
      delete[] cells;      
    }

    bool Valid(const Vectornd& x) const{
      for (int i = 0; i < x.size(); ++i) {
        if (x[i] < xlb[i])
          return false;
        if (x[i] > xub[i])
          return false;
      }
      return true;
    }
   
    int Id(const Vectornd& x) const {

      //      if (n==2) {
      //        return gs[0]*ids[1] + ids[0];
      //      }
      // for n=2
      // id=     gs[0]*ids[1] + ids[0];
      // for n=3
      // id=     gs[0]*gs[1]*ids[2] + gs[0]*ids[1] + ids[0];

      //      assert(Valid(x));
 
      //      int ids[n];  // dimension indices
      int cum = 1; // cumulative offset for next dimension
        
      int id = 0;
      for (int i = 0; i < x.size(); ++i) {
        // index of i-th dimension
        int ind  = floor((x[i] - xlb[i])/(xub[i] - xlb[i])*gs[i]);        
        id += cum*ind;
        if (i < n - 1)
          cum *= gs[i];
      }
      return id;
    } 

    int Index(const Vectornd& x, int i) const {
      return floor((x[i] - xlb[i])/(xub[i] - xlb[i])*gs[i]);        
    }

    Cell<n>* Get(const Vectornd &x, bool checkValid = true) const {
      if (checkValid)
        if  (!Valid(x))
          return 0;
          
      int id = Id(x);
      assert(id >= 0);
      if (id >= nc)
        return 0;
      return cells[id];
    }


    Cell<n>* Get(int id) {
      assert(id >= 0);
      if (id >= nc)
        return 0;
      return cells[id];
    }
    
    Vectornd xlb;  ///< state lower bound
    Vectornd xub;  ///< state upper bound
    Vectorni gs;   ///< number of cells per dimension
    Vectornd cs;   ///< cell length size per dimension

    int nc;  ///< total number of cells
    Cell<n>** cells;
    // or map<Cell<n>*> cells;    
  };
}


#endif
