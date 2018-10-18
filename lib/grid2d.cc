#include "grid2d.h"

using namespace dsl;

Grid2d::Grid2d(int width,
               int height,
               const double* map,
               double sx,
               double sy,
               double max_cost)
    : Grid< Eigen::Vector2d, double >(Eigen::Vector2d(0, 0),
              Eigen::Vector2d(sx * width, sy * height),
              Eigen::Vector2i(width, height)) {
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      int id = j * width + i;
      double cost = map[id]; // cell cost = height/occupany/traversability
      assert(cost >= 0);
      // add this as a cell only if cost is less than a given max cost
      // this is useful if max_cost defines map cells that are untreversable, so
      // they shouldn't be added to the list of cells
      if (cost < max_cost) {
        values[id] = new Cell2d(
            id, Eigen::Vector2d((i + 0.5) * sx, (j + 0.5) * sy), cost);
      }
    }
  }
}

Grid2d::~Grid2d() {
  for (int i = 0; i < nc; ++i) {
    delete values[i];
  }
}
