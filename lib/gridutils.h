#ifndef DSL_GRIDUTILS_H
#define DSL_GRIDUTILS_H

#include <sys/time.h>
#include <Eigen/Dense>
#include "grid.h"
#include "cargrid.h"
#include "cargeom.h"
#include <vector>
#include "utilsimg.h"
#include "utils.h"
#include <thread>

namespace dsl {

/**
 *Load a occupancy grid map from an image file
 * @param filename file to load
 * @param cs cell size
 * @return
 */
Map<bool, 2>::Ptr load(const string& filename, const Eigen::Vector2d &cs);

/**
 * Save an occupancy map with a start, goal and path drawn on it
 * @param map The 2d occupancy map
 * @param filename filename to save the image to
 * @param path pointer to path
 */
void save(const dsl::Map<bool, 2> &map, const string& filename, const std::vector<Eigen::Vector2d> *path = 0);

/**
 * Saves the occupancy map as an image file with primitives drawn on it.
 * @param omap occupancy map
 * @param filename filename to save the image to
 * @param prims primitives
 * @param scale scale to increase resolution of output image
 */
void saveMapWithPrims(const dsl::Map<bool, 2>& omap, std::string filename,
                      const std::vector<vector<Vector2d>>& prims,int scale);
/**
 * Save an occupancy map with a start, goal and the path, as a series of rectangle(car) drawn on it.
 * @param cmap
 * @param filename filename to save the image to
 * @param path the path
 * @param geom geometry of the car
 * @param scale scale to increase resolution of output image
 */
void saveMapWithPath(const dsl::Map<bool, 2>& cmap, std::string filename,
                     const std::vector<Vector3d>& path, const CarGeom& geom, int scale);

/**
 * Save all the slices of an 3d occupancy grid as a series of .ppm file, one for each angle
 * @param cmap
 * @param folder
 * @return
 */
bool saveSlices(const dsl::Map<bool, 3> &cmap, string folder);


/**
  * Takes a 2-dimensional occupancy grid and repeats that for all angles
  * @param omap  2-dimensional occupancy grid
  * @param cmap  3-dimensional occupancy grid
  * @return True if cmap could be modified correctly
  */
 bool MakeSE2Map(const Map<bool, 2> &omap, Map<bool, 3> &cmap);

 /**
  * Takes a 2-dim occupancy grid and uses car geometry to create a 3-dim occ grid for all angles
  * @param geom Geometry of the car
  * @param omap 2-dimensional occupancy grid
  * @param cmap 3-dimensional occupancy grid
  * @return True if cmap could be modified correctly
  */
 bool MakeSE2Map(const CarGeom& geom, const Map<bool, 2> &omap, Map<bool, 3> &cmap, int nthreads=1);

 /**
  * Uses geometry of car to dilate the occupancy grid for a given angle
  * @param dilated the dilated data returned
  * @param omap occupancy map
  * @param geom geometry of car
  * @param theta angle of car
  */
 void DilateMap(vector<bool>& dilated, const Map<bool,2>& omap, const CarGeom& geom, double theta);

// /**
//  * Takes points that make up a kernel and dilates the occupancy grid
//  * @param dmap The dilated map
//  * @param omap The occupancy grid
//  * @param corners Set of points that indicates the bounding rectange of a car
//  */
// void DilateMap(Map<bool,2>& dmap, const Map<bool,2> omap, const vector<Vector2d>& vertices);

}

#endif /* DSL_GRIDUTILS_H */
