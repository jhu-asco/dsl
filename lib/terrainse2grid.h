// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_TERRAINSE2GRID_H_
#define DSL_LIB_TERRAINSE2GRID_H_

#include "grid.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace dsl {
/**
 * Struct to store properties of a cell on terrain. Height, slope(roll pitch), traversibility.
 * The notion of traversibility is correlated to how easy is it to move on a surface. Grass is probably harder
 * to move than say concrete etc.
 */
struct TerrainData{
  using Ptr = std::shared_ptr<TerrainData>;
  TerrainData(double height=std::numeric_limits<double>::quiet_NaN(),
              double traversibility=std::numeric_limits<double>::quiet_NaN());

  double height; ///< Height of the cell relative to the origin.
  double traversibility;///<Traversibility is nan if cell is occupied. Distance*Traversibility is friction loss.
  //double pitch; ///< pitch
  //double roll;  ///< roll

  bool SerializeToString(std::string* str) const;
  bool ParseFromString(const std::string& str);
};

// a cell that stores terrain data along with axy representaion of SE2 pose
using SE2TerrainCell = Cell< Eigen::Vector3d, TerrainData >;
using TerrainSE2GridBase = GridCore< SE2TerrainCell::PointType, SE2TerrainCell::Ptr >;


/**
 * A 3d grid for simple car models with coordinates (theta,x,y). Each cell stores the TerrainData i.e. the height
 * and traversibility of a cell.
 *
 * Authors: Subhransu Mishra
 */
class TerrainSE2Grid : public TerrainSE2GridBase {
public:

  /**
   * Construct terrainse2grid
   * @param cmap Configuration space(angle, x and y) occupancy map
   * @param tmap Terrain data map
   * @param cs cell size of the grid created
   */
  TerrainSE2Grid(const Map<bool, 3> &cmap, const Map<TerrainData,2>& tmap, const Eigen::Vector3d& cs);
};
}

#endif //TERRAINSE2GRID
