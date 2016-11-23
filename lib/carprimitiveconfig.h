// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_CARPRIMITIVECONFIG_H
#define DSL_CARPRIMITIVECONFIG_H
#include <stdint.h>
namespace dsl{
/**
 * Car primitive configuration
 */
struct CarPrimitiveConfig {

  CarPrimitiveConfig(bool fwdonly=true, double tphioverlmx=0.54, double lmin=0.32, double lmax=2,
                  int nl=2, double amax=1.58, int na=12, bool pert=true, bool tocenter=true);

  bool    fwdonly;      //! Decides if the car moves only in the forward direction
  double  tphioverlmax; //! Max(tan(phi)/l) possible for the car
  double  lmin;         //! Min length of the pimitive
  double  lmax;         //! Max length of the primitive
  int     nl;           //! Number of different primitive lengths from lenmin to lenmax. If 1 lmin is chosen.
  double  amax;         //! Maximum angle turned
  int     na;           //! number of steering angles from 0 to phi. If even changed to next odd number
  bool    pert;         //! preturb primitive length or not
  bool    tocenter;     //! moves the primitives to the center of each cell(only in x and y though), hence creates different set of primitives for each angle
};

}

#endif /* DSL_CARPRIMITIVECONFIG_H */
