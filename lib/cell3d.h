#ifndef DSL_CELL3D_H
#define DSL_CELL3D_H

#include <cmath>

namespace dsl {
  class Cell3d {
  public:
    Cell3d(int x = 0, int y = 0, int z = 0, double cost = 0) { 
      p[0] = x;
      p[1] = y;
      p[2] = z;
      this->cost = cost;
    }

    Cell3d(const Cell3d &cell) {
      p[0] = cell.p[0];
      p[1] = cell.p[1];
      p[2] = cell.p[2];
      cost = cell.cost;
    }
   
    double Distance(const Cell3d &cell) const {
      double dx = cell.p[0] - p[0];
      double dy = cell.p[1] - p[1];
      double dz = cell.p[2] - p[2];
      
      return sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    int p[3];
    double cost;
  };
}

#endif
