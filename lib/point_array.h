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
#include <Eigen/Dense>


namespace dsl {

/**
  * An n-dimenensional PointArray consisting of abstract "values"
  * identified by a set of coordinates of type PointT
  * A PointArray provides instant access to the stored values by maintaining
  * an n-dimensional array accessed by convering a PointT to a set of indices.
  * The type ValueT must either be a primitive type or a pointer.
  *
  * Note that this data structure is only viable up to a few dimensions,
  * e.g. dim=5 or 6.
 */
template < typename PointT, typename ValueT >
class PointArray {

  static constexpr int n = PointT::SizeAtCompileTime;

  using Vectorni = Eigen::Matrix< int, n, 1 >;

public:

  PointArray() = delete;

  /**
   * Initialize the array using state lower bound, state upper bound, the number
   * of array values
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of array values per each dimension
   */
  PointArray(const PointT& xlb, const PointT& xub, const Vectorni& gs)
    : xlb(xlb), xub(xub), gs(gs) {
    ds = xub - xlb;
    size = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      mcgs[i] = (i == 0 ? gs[0] : mcgs[i - 1] * gs[i]);
      cs[i] = (xub[i] - xlb[i]) / gs[i];
      gs_div_ds[i] = gs[i] / ds[i];
      size *= gs[i]; // total number of values
    }

    assert(sizeof(ValueT) <= 8);

    values = new ValueT[size];
    memset(values, 0, size * sizeof(ValueT)); // initialize all of them to nil
  }

  /**
   * Initialize the array using state lower bound, state upper bound, the number
   * of array values
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param cs cell dimensions
   */
  PointArray(const PointT& xlb, const PointT& xub, const PointT& cs)
    : xlb(xlb), xub(xub), cs(cs) {
    ds = xub - xlb;
    size = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(cs[i] > 0);
      gs[i] = floor((xub[i] - xlb[i]) / cs[i]);
      assert(gs[i] > 0);
      mcgs[i] = (i == 0 ? gs[0] : mcgs[i - 1] * gs[i]);
      size *= gs[i]; // total number of values
      gs_div_ds[i] = gs[i] / ds[i];
    }

    // here we assume ValueT is of primitive type or is a pointer
    assert(sizeof(ValueT) <= 8);

    values = new ValueT[size];
    memset(values, 0, size * sizeof(ValueT)); // initialize all of them to nil
  }

  PointArray(const PointArray& array)
    : xlb(array.xlb),
      xub(array.xub),
      ds(array.ds),
      gs(array.gs),
      cs(array.cs),
      mcgs(array.mcgs),
      gs_div_ds(array.gs_div_ds),
      size(array.size) {
    values = new ValueT[size];
    memcpy(values, array.values, size * sizeof(ValueT));
  }

  virtual ~PointArray() {
    delete[] values;
  }

  /**
   * Check if point x is whitin array bounds
   * @param x point
   * @return true if within bounds
   */
  virtual bool valid(const PointT& x, double eps = 1e-10) const {
    for (int i = 0; i < x.size(); ++i) {
      if (x[i] < xlb[i] + eps)
        return false;
      if (x[i] > xub[i] - eps)
        return false;
    }
    return true;
  }

  int computeId(const PointT& x) const {
    // unroll loop for n=1,2,3,4 for efficiency
    if (n == 1)
      return index(x, 0);

    if (n == 2)
      return index(x, 0) + mcgs[0] * index(x, 1);

    if (n == 3)
      return index(x, 0) + mcgs[0] * index(x, 1) + mcgs[1] * index(x, 2);

    if (n == 4) {
      return index(x, 0) + mcgs[0] * index(x, 1) + mcgs[1] * index(x, 2) +
        mcgs[2] * index(x, 3);
    }

    int id = 0;
    for (int i = 0; i < x.size(); ++i) {
      id += (i == 0) ? index(x, i) : mcgs[i - 1] * index(x, i);
    }
    return id;
  }


  /**
   * Get the array index of i-th dimension of point x
   * @param x point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(const PointT& x, int i) const {
    return floor((x[i] - xlb[i]) * gs_div_ds[i]);
  }

  /**
   * Get the array index of i-th dimension of point x
   * @param xi i-th coordinate of point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(double xi, int i) const {
    return floor((xi - xlb[i]) * gs_div_ds[i]);
  }


  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return contents of cell
   */
  ValueT data(const PointT& x, bool checkValid = true) const {
    if (checkValid)
      if (!valid(x))
        return empty;

    int id = computeId(x);
    assert(id >= 0);
    if (id >= size)
      return empty;
    return values[id];
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return contents of cell
   */
  void setData(const PointT& x, const ValueT& data, bool checkValid = true) {
    if (checkValid)
      if (!valid(x))
        return;

    int id = computeId(x);
    if (id < 0 || id >= size)
      std::cout << "id=" << id << " x=" << x.transpose()
                << " xlb=" << xlb.transpose() << " xub=" << xub.transpose()
                << std::endl;
    assert(id >= 0 && id < size);
    values[id] = data;
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  ValueT data(int id) const {
    assert(id >= 0);
    if (id >= size)
      return 0;
    return values[id];
  }

  /**
   * Get the cell at a given cell id
   * @param id a non-negative id
   * @return pointer to a cell or 0 if cell does not exist
   */
  void setData(int id, const ValueT& data) {
    assert(id >= 0);
    if (id >= size)
      return;
    values[id] = data;
  }

  static void save(const PointArray< PointT, ValueT >& array, const char* filename) {
    std::ofstream fs(filename, std::fstream::out | std::ios::binary);
    assert(fs.is_open());

    auto index_size = sizeof(xlb[0]);

    for (int i = 0; i < array.xlb.size(); ++i)
      fs.write((char*)&array.xlb[i], index_size);
    for (int i = 0; i < array.xub.size(); ++i)
      fs.write((char*)&array.xub[i], index_size);
    for (int i = 0; i < array.gs.size(); ++i)
      fs.write((char*)&array.gs[i], index_size);

    fs.write((char*)array.values, array.size * sizeof(ValueT));
    fs.close();
  }

  static PointArray< PointT, ValueT >* load(const char* filename) {
    std::ifstream fs(filename, std::fstream::in | std::ios::binary);
    assert(fs.is_open());
    PointT xlb, xub;
    Vectorni gs;

    auto index_size = sizeof(xlb[0]);

    for (int i = 0; i < xlb.size(); ++i)
      fs.read((char*)&xlb[i], index_size);
    for (int i = 0; i < xub.size(); ++i)
      fs.read((char*)&xub[i], index_size);
    for (int i = 0; i < gs.size(); ++i)
      fs.read((char*)&gs[i], index_size);

    PointArray< PointT, ValueT >* array = new PointArray(xlb, xub, gs);
    fs.read((char*)array->values, array->size * sizeof(ValueT));
    fs.close();
    return array;
  }

  PointT xlb; ///< state lower bound
  PointT xub; ///< state upper bound
  PointT ds;  ///< dimensions (ds=xub-xlb)
  Vectorni gs;      ///< number of values per dimension
  PointT cs;  ///< cell length size per dimension
  PointT mcgs;      ///< multiplicative cumulative gs
  PointT gs_div_ds; ///< derived quantity gs/ds

  int size = 0;             ///< total maximum number of values
  ValueT* values = nullptr; ///< array of values

  ValueT empty = 0; ///< empty data
};

}
