#include "grid2d.h"

using namespace dsl;
using namespace Eigen;

Grid2d::Grid2d(int width,
               int height,
               const double* map,
               double sx,
               double sy,
               double maxCost)
    : Grid2dBase(Eigen::Vector2d(0, 0),
              Eigen::Vector2d(sx * width, sy * height),
              Eigen::Vector2i(width, height)) {
  //Iterate over all cells
  auto fun = [&](int id, const Vector2i& gidx){
    Vector2d cc; //CellCenter
    bool gotcenter = CellCenter(gidx, &cc);
    assert(gotcenter);
    double cost = map[id]; // cell cost = height/occupany/traversability
    assert(cost >= 0);
    if (cost < maxCost) // add this as a cell only if cost is less than a given max cost
      this->set_cells(id,XyCostCell(id, cc, cost));
  };
  LoopOver(fun);
}
