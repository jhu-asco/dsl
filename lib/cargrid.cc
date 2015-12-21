#include "cargrid.h"

using namespace dsl;
using namespace Eigen;

CarGrid::CarGrid(int width, int height, double *map, 
                 double sx, double sy, double sa, double costScale,
                 double maxCost) :
  Grid<3>(Vector3d(0,0,-M_PI + sa/2), Vector3d(sx*width, sy*height, M_PI + sa/2), 
          Vector3i(width, height, (int)round(2*M_PI/sa))), maxCost(maxCost) {  

  const int &angRes = gs[2];

  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {      
      int mid = j*width + i;
      double cost = map[mid]*costScale; // cell cost = height/occupany/traversability      
      for (int k = 0; k < angRes; ++k) {      
        int id = k*width*height + j*width + i;
        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so 
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[id] = new Cell<3>(xlb + Vector3d((i + 0.5)*sx, (j + 0.5)*sy, (k + 0.5)*sa), 
                                  Vector3d(sx/2, sy/2, sa/2), cost);
        }
      }
    }
  }
}
