// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_LIB_PPM_READER_H_
#define DSL_LIB_PPM_READER_H_
#include <string>
#include <fstream>
#include <vector>

namespace dsl{
/**
 * Container that can store RGB image with a maximum bit depth of 0xffff
 */
struct ImageRGB{
  enum BitDepth: uint16_t{ BD8  = 0xff, BD16 = 0xffff };

  int w, h;

  uint16_t bitdepth; //supported values: 1, 0xff, 0xffff
  std::vector<uint16_t> rdata;
  std::vector<uint16_t> gdata;
  std::vector<uint16_t> bdata;

  void ChangeBitDepth(BitDepth bitdepth_new);
};

/**
 * Loads a .ppm file to an RGB image
 * @param img The iamge retrieved
 * @param filename Filename with path
 * @return True if image could be loaded
 */
bool LoadPpm(ImageRGB& img, const std::string& filename);

/**
 * Save a RGB image to .ppm file
 * @param img
 * @param filename
 * @return
 */
bool SavePpm(ImageRGB& img,  const std::string& filename);

/**
 * Removes the comment section of a std::ifstream
 * @param f file stream
 */
void RemoveComment(std::ifstream &f);

} //namespace dsl
#endif /* DSL_LIB_PPM_READER_H_ */
