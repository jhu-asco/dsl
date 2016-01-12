// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LINECONNECTIVITY_H
#define DSL_LINECONNECTIVITY_H

#include "gridconnectivity.h"
#include "grid.h"

namespace dsl {
  
  using namespace std;
  
  /**
   * The most basic grid connectivity: each cell is connected with "lines" to
   * other cells, and each line has an associate cost.
   *
   * Author: Marin Kobilarov
   */
  template<int n>
    class LineConnectivity : public GridConnectivity<n> {
  public:
    
    typedef Matrix<double, n, 1> Vectornd;
    
    /**
     * Initialize connectivity using a grid
     * @param grid the grid
     */
    LineConnectivity(const Grid<n> &grid);

    /**
     * Initialize connectivity using a grid, lines, and costs
     * @param grid the grid
     * @param lines the lines
     * @param costs the costs
     */    
    LineConnectivity(const Grid<n> &grid,
                     const vector<Vectornd> &lines,
                     const vector<double> &costs);
    

    virtual bool operator()(const Cell<n>& from, 
                            std::vector<GridPath<n> >& paths, 
                            bool fwd = true) const;

    const Grid<n> &grid;          ///< grid    
    std::vector<Vectornd> lines;  ///< line vectors (directions) connecting to other cells    
    std::vector<double> costs;    ///< cost along each direction
    
  };
  
  
  template<int n>
    LineConnectivity<n>::LineConnectivity(const Grid<n> &grid) : grid(grid) {
  }
  
  template<int n>
    LineConnectivity<n>::LineConnectivity(const Grid<n> &grid, 
                                          const vector<Vectornd> &lines,
                                          const vector<double> &costs) : grid(grid), lines(lines), costs(costs) {
  }
  
  template<int n>
    bool LineConnectivity<n>::operator()(const Cell<n>& from, 
                                         vector<GridPath<n> > &paths,
                                         bool fwd) const {
    
    paths.clear();
    for (int i = 0; i < (int)(lines.size()); ++i) {
      GridPath<n> path;
      path.cells.push_back(from);
      
      Vectornd x = from.c;
      if (fwd)
        x += lines[i];
      else
        x -= lines[i];
      
      if (!grid.Valid(x))
        continue;
      
      const Cell<n>* to = grid.Get(x, false);
      if (!to) // cell might be empty
        continue;
      
      path.cells.push_back(*to);    
      path.cost = (1 + (from.cost + to->cost)/2)*this->costs[i];
      path.fwd = fwd;
      paths.push_back(path);
    }
    return true;
  } 
}

#endif
