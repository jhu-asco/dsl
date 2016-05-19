// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Subhransu Mishra <subhransu.kumar.mishra@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "carprimcfg.h"
using namespace dsl;

CarPrimitiveCfg::CarPrimitiveCfg(bool fwdonly, double tphioverlmx, double lmin, double lmax,
                                 uint16_t nl, double amax, uint16_t na, bool pert, bool tocenter):
        fwdonly(fwdonly), tphioverlmax(tphioverlmx), lmin(lmin),lmax(lmax),
        nl(nl), amax(amax), na(na), pert(pert), tocenter(tocenter){

}

