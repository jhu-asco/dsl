// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CELL_H
#define DSL_CELL_H

#include <Eigen/Dense>

namespace dsl {

  using namespace Eigen;
  
  template<int n>
    class Cell {
    
    typedef Matrix<double, n, 1> Vectornd;
    
  public:

  Cell() : cost(0) {
    }
    
  Cell(const Vectornd &c, double cost = 1) : c(c), cost(cost) {
    }
    
  Cell(const Vectornd &c, const Vectornd &r, double cost = 1) : c(c), r(r), cost(cost) {
    }
    
  Cell(const Cell &cell) : c(cell.c), r(cell.r), cost(cell.cost) {      
    }
    
    
    bool Inside(const Vectornd &x) {
      for (int i = 0; i < x.size(); ++i)
        if (x[i] < c[i] - r[i] ||
            x[i] > c[i] + r[i])
          return false;
      return true;
    }
    
    Vectornd c;   ///< center of cell
    Vectornd r;   ///< half-distance of each cell side
    
    double cost;  ///< cost of cell (typically this is 0 if unoccupied) 
  };

}


#endif
