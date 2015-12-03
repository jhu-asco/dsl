#ifndef DSL_CELL2D_H
#define DSL_CELL2D_H

#include <cmath>

namespace dsl {
  class Cell2d {
  public:
    Cell2d(int x = 0, int y = 0, double cost = 0) { 
      p[0] = x;
      p[1] = y;
      this->cost = cost;
    }

    Cell2d(const Cell2d &cell) {
      p[0] = cell.p[0];
      p[1] = cell.p[1];
      cost = cell.cost;
    }
   
    double Distance(const Cell2d &cell) {
      double dx = cell.p[0] - p[0];
      double dy = cell.p[1] - p[1];
      
      return sqrt(dx*dx + dy*dy);
    }
    
    int p[2];
    double cost;
  };
}

#endif
