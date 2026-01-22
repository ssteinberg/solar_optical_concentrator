#pragma once

#include "linalg.h"
using namespace linalg::aliases;

// floating-point precision
using float_type = double;

// doinking intersection points
constexpr bool DOINKING = true;
constexpr float_type DOINK = 1e-5;

// finite system

#define BUILD_FINITE                                true
#define BUILD_FINITE_ELLIPSE                        true
#define BUILD_FINITE_PARABOLA                       false

// infinite system
#define BUILD_INFINITE                              false
#define BUILD_INFINITE_ARBITRARY                    false

// source and target
constexpr bool FLAT_SOURCE = true;
constexpr bool FLAT_TARGET = false;
constexpr bool IGNORE_SOURCE = true;
