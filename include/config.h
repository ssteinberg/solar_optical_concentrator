#pragma once

#include "shape.h"

#include "linalg.h"
using namespace linalg::aliases;

// floating-point precision
using float_type = double;

// doinking intersection points
constexpr bool DOINKING = true;
constexpr float_type DOINK = 1e-5;

// finite system
#define BUILD_FINITE                                false
#define BUILD_FINITE_ELLIPSE                        false
#define BUILD_FINITE_PARABOLA                       false

// infinite system
#define BUILD_INFINITE                              true
#define BUILD_INFINITE_ARBITRARY                    false

// source and target
constexpr bool FLAT_SOURCE = true;
constexpr bool IGNORE_SOURCE = true;
constexpr Shape TARGET_SHAPE = Shape::FLAT;
constexpr float_type ELLIPTICAL_TARGET_X_RADIUS = 0.001;
constexpr float_type ELLIPTICAL_TARGET_Y_RADIUS = 0.005;

// S(β)
constexpr bool SHOULD_USE_MINIFICATION_CONSTANT = false;
