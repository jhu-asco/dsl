#include "grid3d.h"

using namespace dsl;

Grid3d::Grid3d(int length,
               int width,
               int height,
               const double* map,
               double sx,
               double sy,
               double sz,
               double costScale,
               double maxCost)
  : Grid< 3 >(Vector3d(0, 0, 0),
              Vector3d(sx * length, sy * width, sz * height),
              Vector3i(length, width, height)) {
  for (int i = 0; i < length; ++i) {
    for (int j = 0; j < width; ++j) {
      for (int k = 0; k < height; ++k) {
        int id = k * length * width + j * length + i;
        double cost =
            map[id] * costScale; // cell cost = height/occupany/traversability
        assert(cost >= 0);
        // add this as a cell only if cost is less than a given max cost
        // this is useful if maxCost defines map cells that are untreversable,
        // so
        // they shouldn't be added to the list of cells
        if (cost < maxCost) {
          cells[id] = new Cell< 3 >(
              Eigen::Vector3d((i + 0.5) * sx, (j + 0.5) * sy, (k + 0.5) * sz),
              Eigen::Vector3d(sx / 2, sy / 2, sz / 2),
              cost);
        }
      }
    }
  }
}
