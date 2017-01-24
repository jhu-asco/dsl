#include "grid2d.h"

using namespace dsl;

Grid2d::Grid2d(int width,
               int height,
               const double* map,
               double sx,
               double sy,
               double maxCost)
    : Grid< Eigen::Vector2d, double >(Eigen::Vector2d(0, 0),
              Eigen::Vector2d(sx * width, sy * height),
              Eigen::Vector2i(width, height)) {
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      int id = j * width + i;
      double cost = map[id]; // cell cost = height/occupany/traversability
      assert(cost >= 0);
      // add this as a cell only if cost is less than a given max cost
      // this is useful if maxCost defines map cells that are untreversable, so
      // they shouldn't be added to the list of cells
      if (cost < maxCost) {
        cells_[id].reset(new Cell< Eigen::Vector2d, double>(id,
                                                       Eigen::Vector2d((i + 0.5) * sx, (j + 0.5) * sy),
                                                       cost));
      }
    }
  }
}
