// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Eigen/Dense>

namespace dsl {

/**
 * An n-dimenensional IndexedArray consisting of abstract "values", or elements
 * identified by a set of coordinates of type PointT, each cell
 * containing data of type DataT.
 * A IndexedArray provides instant access to the elements by maintaining
 * an n-dimensional array of pointers to values. values that are empty,
 * e.g. that are inside obstacles, correspond to null pointers.
 *
 * Note that this data structure is only viable up to a few dimensions,
 * e.g. dim=5 or 6.
 */
template < class PointT, class ValueT >
struct IndexedArray {
  using Vectorni = Eigen::Matrix< int, PointT::SizeAtCompileTime, 1 >;

  IndexedArray() = delete;

  /**
   * Initialize the IndexedArray using state lower bound, state upper bound, the
   * number
   * of IndexedArray values
   * @param xlb state lower bound
   * @param xub state upper bound
   * @param gs number of IndexedArray values per each dimension
   */
  IndexedArray(const PointT& xlb, const PointT& xub, const Vectorni& gs)
    : n(xlb.size()), xlb(xlb), xub(xub), gs(gs) {
    ds = xub - xlb;
    size = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] <= xub[i]);
      assert(gs[i] > 0);
      size *= gs[i]; // total number of values
      mcgs[i] = (i == 0 ? gs[0] : mcgs[i - 1] * gs[i]);
      cs[i] = (xub[i] - xlb[i]) / gs[i];
      gs_div_ds[i] = 1.0 / cs[i];
    }
    values = new ValueT[size];

    // assume ValueT is of primitive type or is a pointer
    assert(sizeof(ValueT) <= 8);
    memset(values, 0, size * sizeof(ValueT)); // initialize all of them nil
  }

  /**
 * Initialize the map using state lower bound, state upper bound, the number
 * of map values
 * @param xlb state lower bound
 * @param xub state upper bound
 * @param cs cell dimensions
 */
  IndexedArray(const PointT& xlb, const PointT& xub, const PointT& cs)
    : n(xlb.size()), xlb(xlb), xub(xub), cs(cs) {
    ds = xub - xlb;
    size = 1;
    for (int i = 0; i < n; ++i) {
      assert(xlb[i] < xub[i]);
      assert(cs[i] > 0);
      gs[i] = floor((xub[i] - xlb[i]) / cs[i]);
      mcgs[i] = (i == 0 ? gs[0] : mcgs[i - 1] * gs[i]);
      size *= gs[i]; // total number of values
      gs_div_ds[i] = gs[i] / ds[i];
    }
    values = new ValueT[size];

    // here we assume ValueT is of primitive type or is a pointer
    assert(sizeof(ValueT) <= 8);
    memset(values, 0, size * sizeof(ValueT)); // initialize all of them nil
  }

  IndexedArray(const IndexedArray& array)
    : n(array.n),
      xlb(array.xlb),
      xub(array.xub),
      ds(array.ds),
      cs(array.cs),
      gs(array.gs),
      mcgs(array.mcgs),
      gs_div_ds(array.gs_div_ds),
      size(array.size) {
    values = new ValueT[size];
    memcpy(values, array.values, size * sizeof(ValueT));
  }

  virtual ~IndexedArray() {
    delete[] values;
  }

  /**
   * Check if point x is within IndexedArray bounds
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

  /**
   * Get an id of point x useful for direct lookup in the IndexedArray array
   * @param x point
   * @return a computed id
   */
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
      id += mcgs[i] * index(x, i);
    }
    return id;
  }

  /**
   * Get the IndexedArray index of i-th dimension of point x
   * @param x point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(const PointT& x, int i) const {
    return floor((x[i] - xlb[i]) / gs_div_ds[i]);
  }

  /**
   * Get the IndexedArray index of i-th dimension of point x
   * @param xi i-th coordinate of point
   * @param i coordinate index
   * @return index into cell array
   */
  int index(double xi, int i) const {
    return floor((xi - xlb[i]) / gs_div_ds[i]);
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
    // if (id<0 || id>=size)
    //  std::cout << "id=" << id << " x=" << x.transpose() << " xlb=" <<
    //  xlb.transpose() << " xub=" << xub.transpose() << std::endl;
    assert(id >= 0 && id < size);
    values[id] = data;
  }

  /**
   * Get the cell at position x
   * @param x point
   * @param checkValid whether to check if within valid bounds (more efficient
   * if checkValid=0 but dangerous)
   * @return pointer to a cell or 0 if cell does not exist
   */
  ValueT data(const PointT& x, bool checkValid = true) const {
    if (checkValid)
      if (!valid(x))
        return 0;
    return data(computeId(x));
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

  int n; ///< IndexedArray dimension

  PointT xlb;  ///< state lower bound
  PointT xub;  ///< state upper bound
  PointT ds;   ///< dimensions (ds=xub-xlb)
  PointT cs;   ///< cell length size per dimension
  Vectorni gs; ///< number of values per dimension

  PointT mcgs;      ///< multiplicative cumulative gs
  PointT gs_div_ds; ///< derived quantity gs/ds

  int size = 0;             ///< number of values in IndexedArray
  ValueT* values = nullptr; ///< array of data
};
}
