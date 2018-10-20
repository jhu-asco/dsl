// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <fstream>

namespace dsl {

/**
 * An n-dimenensional occupancy map storing type T in every cell
 */
template < typename T, int n >
class Map {
  using Vectornd = Eigen::Matrix< double, n, 1 >;
  using Vectorni = Eigen::Matrix< int, n, 1 >;

public:
  /**
   * Initialize the map using state lower bound, state upper bound, the number
   * of map cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of map cells per each dimension
   */
  Map(const Vectornd& xlb, const Vectornd& xub, const Vectorni& gs)
    : xlb(xlb), xub(xub), gs(gs) {
    ds = xub - xlb;
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      nc *= gs[i]; // total number of cells
      cs[i] = (xub[i] - xlb[i]) / gs[i];
    }

    cells = new T[nc];
    memset(cells, 0, nc * sizeof(T)); // initialize all of them to nil
  }

  /**
   * Initialize the map using state lower bound, state upper bound, the number
   * of map cells
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param cs cell dimensions
   */
  Map(const Vectornd& xlb, const Vectornd& xub, const Vectornd& cs)
    : xlb(xlb), xub(xub), cs(cs) {
    ds = xub - xlb;
    nc = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(cs[i] > 0);
      gs[i] = floor((xub[i] - xlb[i]) / cs[i]);
      nc *= gs[i]; // total number of cells
    }

    cells = new T[nc];
    memset(cells, 0, nc * sizeof(T)); // initialize all of them to nil
  }

  Map(const Map& map)
    : xlb(map.xlb),
      xub(map.xub),
      ds(map.ds),
      gs(map.gs),
      cs(map.cs),
      nc(map.nc) {
    cells = new T[nc];
    memcpy(cells, map.cells, nc * sizeof(T));
  }

  virtual ~Map() {
    delete[] cells;
  }

  /**
   * Check if point x is whitin map bounds
   * @param x point
   * @return true if within bounds
   */
  virtual bool valid(const Vectornd& x, double eps = 1e-10) const {
    for (int i = 0; i < x.size(); ++i) {
      if (x[i] < xlb[i] + eps)
        return false;
      if (x[i] > xub[i] - eps)
        return false;
    }
    return true;
  }

  /**
   * Get an id of point x useful for direct lookup in the map array
   * @param x point
   * @return a computed id
   */
  int computeId(const Vectornd& x) const {
    // unroll loop for n=1,2,3,4 for efficiency
    if (n == 1)
      return index(x, 0);

    if (n == 2)
      return index(x, 0) + gs[0] * index(x, 1);

    if (n == 3)
      return index(x, 0) + gs[0] * index(x, 1) + gs[0] * gs[1] * index(x, 2);

    if (n == 4) {
      int cum = gs[0] * gs[1];
      return index(x, 0) + gs[0] * index(x, 1) + cum * index(x, 2) +
          cum * gs[2] * index(x, 3);
    }

    // general for any n
    int cum = 1; // cumulative offset for next dimension

    int id = 0;
    for (int i = 0; i < x.size(); ++i) {
      // index of i-th dimension
      id += cum * index(x, i);
      if (i < n - 1)
        cum *= gs[i];
    }
    return id;
  }

  /**
   * Get the map index of i-th dimension of point x
   * @param x point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(const Vectornd& x, int i) const {
    return floor((x[i] - xlb[i]) / ds[i] * gs[i]);
  }

  /**
   * Get the map index of i-th dimension of point x
   * @param xi i-th coordinate of point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(double xi, int i) const {
    return floor((xi - xlb[i]) / ds[i] * gs[i]);
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return contents of cell
   */
  T data(const Vectornd& x, bool checkValid = true) const {
    if (checkValid)
      if (!valid(x))
        return empty;

    int id = computeId(x);
    assert(id >= 0);
    if (id >= nc)
      return empty;
    return cells[id];
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return contents of cell
   */
  void setData(const Vectornd& x, const T& data, bool checkValid = true) {
    if (checkValid)
      if (!valid(x))
        return;

    int id = computeId(x);
    if (id < 0 || id >= nc)
      std::cout << "id=" << id << " x=" << x.transpose()
                << " xlb=" << xlb.transpose() << " xub=" << xub.transpose()
                << std::endl;
    assert(id >= 0 && id < nc);
    cells[id] = data;
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  T data(int id) const {
    assert(id >= 0);
    if (id >= nc)
      return 0;
    return cells[id];
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  void setData(int id, const T& data) {
    assert(id >= 0);
    if (id >= nc)
      return;
    cells[id] = data;
  }

  static void save(const Map< T, n >& map, const char* filename) {
    std::ofstream fs(filename, std::fstream::out | std::ios::binary);
    assert(fs.is_open());

    for (int i = 0; i < map.xlb.size(); ++i)
      fs.write((char*)&map.xlb[i], sizeof(double));
    for (int i = 0; i < map.xub.size(); ++i)
      fs.write((char*)&map.xub[i], sizeof(double));
    for (int i = 0; i < map.gs.size(); ++i)
      fs.write((char*)&map.gs[i], sizeof(double));

    fs.write((char*)map.cells, map.nc * sizeof(T));
    fs.close();
  }

  static Map< T, n >* load(const char* filename) {
    std::ifstream fs(filename, std::fstream::in | std::ios::binary);
    assert(fs.is_open());
    Vectornd xlb, xub;
    Vectorni gs;

    for (int i = 0; i < xlb.size(); ++i)
      fs.read((char*)&xlb[i], sizeof(double));
    for (int i = 0; i < xub.size(); ++i)
      fs.read((char*)&xub[i], sizeof(double));
    for (int i = 0; i < gs.size(); ++i)
      fs.read((char*)&gs[i], sizeof(double));

    Map< T, n >* map = new Map(xlb, xub, gs);
    fs.read((char*)map->cells, map->nc * sizeof(T));
    fs.close();
    return map;
  }

  Vectornd xlb; ///< state lower bound
  Vectornd xub; ///< state upper bound
  Vectornd ds;  ///< dimensions (ds=xub-xlb)
  Vectorni gs;  ///< number of cells per dimension
  Vectornd cs;  ///< cell length size per dimension

  int nc = 0;         ///< total maximum number of cells
  T* cells = nullptr; ///< array of cells

  T empty = 0; ///< empty data
};
}
