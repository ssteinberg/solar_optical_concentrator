#pragma once

#include "types.h"

struct HitInfo {
    // regular stuff
    float_type l;
    float_vec p, n;
    Type t;

    // length up mirror
    float_type ml;

    // arbitrary target
    float_vec u;
    float_type rtraced, rmin, rmax, rmean;
    float_type rmin2, rmax2, rmean2;
    float_type a, sina, tana;
};
