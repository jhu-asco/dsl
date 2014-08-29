// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_COST_H
#define DSL_COST_H

#include "vertex.h"


namespace dsl {
  /**
   * Interface defining the cost b/n two vertices.
   * Graph search algorithms often require two types of cost:
   *
   * * the best possible real cost that a path b/n two vertices can have
   *
   * * a heuristic cost that is an underestimate of that real cost but
   *     is useful during path searching.
   *
   * This interface defines the two types of costs. The heuristic function
   * is more important and subclasses of this interface are required to
   * provide it. The real function is used only for tie-breaking and
   * is not crucial but could change performance in certain cases.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  class Cost {
  public:
    /**
     * Heuristic distance 
     * (should be admissible, i.e. less than the real minimum distance)
     * b/n two vertices.
     * Subclasses must provide this function
     * @param va first vertex
     * @param vb second vertex
     * @return heuristic distance (optimal cost)
     */
    virtual double Heur(const Vertex &va, const Vertex &vb) const = 0;
    
    
    /**
     * Real best possible distance b/n two vertices.
     * Subclasses should optionally provide this function.
     * By default it is the heuristic function + epsilon 
     * @param va first vertex
     * @param vb second vertex
     * @return real minimum possible cost b/n va and vb
     */
    virtual double Real(const Vertex &va, const Vertex &vb) const;
  };
}

#endif
