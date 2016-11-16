#ifndef DSL_GRIDUTILS_H
#define DSL_GRIDUTILS_H

#include <sys/time.h>
#include <Eigen/Dense>
#include "grid.h"
#include "cargrid.h"
#include <vector>
#include "utilsimg.h"
#include "utils.h"

namespace dsl {

Map<bool, 2> load(const char* filename, const Eigen::Vector2d &cs);

/**
 * Save an occupancy with a start, goal and path drawn on it
 * @param map
 * @param filename
 * @param path
 * @param start
 * @param goal
 */
void save(const dsl::Map<bool, 2> &map, const char* filename, const std::vector<Eigen::Vector2d> *path = 0);

void saveMapWithPath(const dsl::Map<bool, 2>& cmap, std::string filename,
                     const std::vector<Vector3d>& path, const CarGeom& geom, int scale);


/**
 * Get the corners of car which is rotate by an angle theta. The geometry of car is defined by geom.
 * Corners in the order (right-back, right-front, left-front and left-back)
 * @param verts2d_rotd_pix
 * @param geom
 * @param theta
 */
void getCarCorners(Matrix2x4d& verts2d_rotd_pix,const CarGeom& geom, double theta);

}

#endif /* DSL_GRIDUTILS_H */
