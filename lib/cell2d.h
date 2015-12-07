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

    /**
    * Operator to access points
    * @param index : Index of the elemen accessing 0 or 1
    * @return returns the the value of the cell at the index. If out of bounds, returns -1.
    */
    int operator [](int index)
    {
      if(index == 0 || index == 1)
      {
        return p[index];
      }
      return -1;
    }
    
    int p[2];
    double cost;
  };
}

#endif
