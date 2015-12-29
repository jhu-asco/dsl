#include "cargrid.h"
#include "utils.h"

using namespace dsl;
using namespace Eigen;

CarGrid::CarGrid(int width, int height, double *map, 
                 double sx, double sy, double sa, double costScale,
                 double maxCost) :
  Grid<3, Matrix3d>(Vector3d(-M_PI + sa/2, 0, 0), Vector3d(M_PI + sa/2, sx*width, sy*height), 
                    Vector3i((int)round(2*M_PI/sa), width, height)), maxCost(maxCost) {  
  
  const int &angRes = gs[0];
  /*
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
          cells[id] = new Cell<3>(xlb + Vector3d((k + 0.5)*sa, (i + 0.5)*sx, (j + 0.5)*sy), 
                                  Vector3d(sa/2, sx/2, sy/2), cost);
        }
      }
    }
  }
  */


  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {      
      int mid = j*width + i;
      double cost = map[mid]*costScale; // cell cost = height/occupany/traversability      
      for (int k = 0; k < angRes; ++k) {      
        int id = j*angRes*width + i*angRes + k;

        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable, so 
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[id] = new SE2Cell(xlb + Vector3d((k + 0.5)*sa, (i + 0.5)*sx, (j + 0.5)*sy), 
                                  Vector3d(sa/2, sx/2, sy/2), cost);
          se2_q2g(cells[id]->data, cells[id]->c);
        }
      }
    }
  }

}
