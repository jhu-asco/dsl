#include "grid3d.h"

using namespace dsl;
using namespace Eigen;

Grid3d::Grid3d(int length,
               int width,
               int height,
               const double* map,
               double sx,
               double sy,
               double sz,
               double costScale,
               double maxCost)
    : Grid3dBase(Eigen::Vector3d(0, 0, 0),
      Eigen::Vector3d(sx * length, sy * width, sz * height),
      Eigen::Vector3i(length, width, height)) {
  //Iterate over all cells
  auto fun = [&](int id, const Vector3i& gidx){
    Vector3d cc; //CellCenter
    bool gotcenter = CellCenter(gidx, &cc);
    assert(gotcenter);
    double cost = map[id] * costScale; // cell cost = height/occupany/traversability
    assert(cost >= 0);
//    if (cost < maxCost) // add this as a cell only if cost is less than a given max cost
//      this->set_cells(id,XyzCostCell(id, cc, cost));
  };
  LoopOver(fun);
}
