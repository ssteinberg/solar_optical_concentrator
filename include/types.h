#pragma once

#include "config.h"

typedef linalg::vec<float_type, 2> float_vec;

constexpr float_type PI = 3.14159265358979;
constexpr float_type PI_OVER_TWO = PI / 2;
constexpr float_type ONE_OVER_PI = 1 / PI;
constexpr float_type DEG_TO_RAD = PI / 180;
constexpr float_type RAD_TO_DEG = 180 / PI;
constexpr float_type EPSILON = 0.05;
constexpr float_type EPSILON_OVER_TWO = EPSILON / 2;

enum struct Shape {
    FLAT,
    CYLINDRICAL,
    ELLIPTICAL,
    PARABOLIC,
    CONSTRUCTED
};

enum struct Type {
    SOURCE,
    TARGET,
    MIRROR,
    MIRROR1,
    MIRROR2,
    MIRROR_1a,
    MIRROR_1b,
    MIRROR_2a,
    MIRROR_2b,
    BARRIER
};
