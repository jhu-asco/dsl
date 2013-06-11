#ifndef DSL_GRIDCOST_H
#define DSL_GRIDCOST_H

#include "cost.h"

namespace dsl {
  /**
   * Grid cost interface.
   *
   * Author: Marin Kobilarov -- Copyright (C) 2004
   */
  class GridCost : public Cost {
  public:
    double Heur(const Vertex &va, const Vertex &vb) const;       
    double Real(const Vertex &va, const Vertex &vb) const;    
  };
}

#endif
