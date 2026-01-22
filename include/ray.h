#pragma once

#include "types.h"

struct Ray {
    float_vec o, d;
    Ray() : o(), d(float_vec(0.0, 1.0)) {}
    Ray(const float_vec& o, const float_vec& d) : o(o), d(d) {}
};
