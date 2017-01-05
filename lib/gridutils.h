#ifndef DSL_GRIDUTILS_H
#define DSL_GRIDUTILS_H

#include "grid.h"
#include "terrainse2grid.h"
#include "cargeom.h"
#include <vector>
#include <sys/time.h>
#include <Eigen/Dense>

namespace dsl {
/**
 * Load a 2D occupancy grid map from an image file(.ppm)
 * @param filename file to load
 * @param cs cell size
 * @return The loaded map. Nullptr if unsuccesful.
 */
Map<bool, 2>::Ptr LoadPpm(const std::string& filename, const Vector2d &cs);

/**
 * Save a 2D occupancy map to a .ppm file
 * @param map The 2D occupancy map
 * @param filename filename to save the image to
 */
bool SavePpm(const dsl::Map<bool, 2> &map, const std::string& filename);

/**
 * Save all the slices of an 3D occupancy grid as a series of .ppm file, one for each angle
 * @param cmap Occupancy grid in angle x and y
 * @param folder folder in which to put the images
 * @return True if saving was successful.
 */
bool SavePpm(const dsl::Map<bool, 3> &cmap, std::string folder);

/**
 * Load an occupancy map from .omap file
 * @param omapfile filename with .omap extension
 * @return The loaded map. Nullptr if unsuccesful.
 */
Map<bool, 2>::Ptr LoadOmap(const std::string& omapfile);

/**
 * Load a terrain map from .tmap file
 * @param tmapfile filename with .tmap extension
 * @return The loaded map. Nullptr if unsuccesful.
 */
Map<TerrainData, 2>::Ptr LoadTmap(const std::string& tmapfile);

///**
// * Save an occupancy map to .omap file
// * @param omap The occupancy map to be saved
// * @param omapfile filename with .omap extension
// * @return True if saving was successful
// */
//bool saveMap(Map<bool, 2>::Ptr& omap, const string& omapfile);
//
///**
// * Save a terrain map to .tmap file
// * @param tmap The terrain map to be saved
// * @param tmapfile filename with .tmap extension
// * @return True if saving was successful
// */
//bool saveMap(Map<TerrainData, 2>::Ptr& tmap, const string& tmapfile);

/**
 * Save an occupancy map as .ppm image with a start(green), goal(red) and the waypoints(blue) as points.
 * @param omap 2D occupancy map
 * @param filename .ppm filename to save the image to
 * @param path 2D path(x and y)
 */
void SavePpmWithPath(const dsl::Map<bool, 2> &map, const std::string& filename, const std::vector<Vector2d>& path);

/**
 * Save an occupancy map as .ppm image with a start(green), goal(red) and the waypoints(blue) as points.
 * If geometry is available then it plots rectangles instead.
 * @param omap 2D occupancy map
 * @param filename .ppm filename to save the image to
 * @param scale scale to increase resolution of output image for clarity.
 * @param path 3D path(angle, x and y)
 * @param geom pointer to optional geometry of the car
 */
bool SavePpmWithPath(const dsl::Map<bool, 2>& omap, std::string filename, int scale,
                      const std::vector<Vector3d>& path, const CarGeom* geom = 0);

/**
 * Save an terrain map as .ppm image with a start(green), goal(red) and the waypoints(blue) as points.
 * If geometry is available then it plots rectangles instead.
 * @param tmap 2D map with terrain data
 * @param filename .ppm filename to save the image to
 * @param scale scale to increase resolution of output image for clarity.
 * @param path path
 * @param geom pointer to optional geometry of the car
 */
bool SavePpmWithPath(const dsl::Map<TerrainData, 2>& tmap, std::string filename, int scale,
                     const std::vector<Vector3d>& path, const CarGeom* geom = 0);

/**
 * Save an occupancy map as .ppm image with a set of motion primitives at start location
 * @param omap 2D occupancy map
 * @param filename .ppm filename to save the image to
 * @param scale scale to increase resolution of output image for clarity.
 * @param prims primitives
 */
bool SavePpmWithPrimitives(const dsl::Map<bool, 2>& omap, std::string filename, int scale,
                      const std::vector<vector<Vector2d>>& prims);

/**
 * Saves the terrain map as .ppm image with a set of motion primitives at start location
 * @param tmap terrain map
 * @param filename filename to save the image to
 * @param scale scale to increase resolution of output image for clarity
 * @param prims primitives
 */
bool SavePpmWithPrimitives(const dsl::Map<TerrainData, 2>& tmap, std::string filename,int scale,
                      const std::vector<vector<Vector2d>>& prims);

/**
 * Takes a 2-dimensional occupancy map and repeats that for all angles to create configuration
  * space(a,x,y) occupancy map
 * @param omap  2-dimensional occupancy map
 * @param csa cell size for angle
 * @return configuration space(yaw, x and y) occupancy map
 */
 Map<bool, 3>::Ptr MakeCmap(const Map<bool, 2> &omap, double csa);

 /**
  * Takes a 2-dim occupancy grid and car geometry to create a configuration space(a,x,y)
  * @param omap 2-dimensional map of occupancy data
  * @param csa cell size for angle
  * @param geom Geometry of the car
  * @param nthreads number of threads for making cmap
  * @return configuration space(yaw, x and y) occupancy map
  */
 Map<bool, 3>::Ptr MakeCmap(const Map<bool, 2> &omap,double csa, const CarGeom& geom, int nthreads = 1);

 /**
  * Takes a map(in x and y) of terrain data and creates a configuration space occupancy map
  * @param tmap 2-dimensional map of terrain data
  * @param csa cell size for angle
  * @return configuration space(yaw, x and y) occupancy map
  */
  Map<bool, 3>::Ptr MakeCmap(const Map<TerrainData, 2> &tmap, double csa);

  /**
   * Takes a map(in x and y) of terrain data and car geometry to create a configuration
   * space occupancy map
   * @param tmap 2-dimensional map of terrain data
   * @param csa cell size for angle
   * @param geom geometry of car
   * @param nthreads number of threads for making cmap
   * @return configuration space(yaw, x and y) occupancy map
   */
  Map<bool, 3>::Ptr MakeCmap(const Map<TerrainData, 2> &tmap,double csa, const CarGeom& geom, int nthreads = 1);

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
