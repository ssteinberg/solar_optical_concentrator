#pragma once

#include "linalg.h"
using namespace linalg::aliases;

// floating-point precision
using float_type = double;

// doinking intersection points
constexpr bool DOINKING = true;
constexpr float_type DOINK = 1e-5;

// finite system
constexpr bool BUILD_FINITE = true;
constexpr bool BUILD_FINITE_ELLIPSE = true;
constexpr bool BUILD_FINITE_PARABOLA = false;

// infinite system
constexpr bool BUILD_INFINITE = false;
constexpr bool BUILD_INFINITE_ARBITRARY = false;

// source and target
constexpr bool FLAT_SOURCE = true;
constexpr bool FLAT_TARGET = false;
constexpr bool IGNORE_SOURCE = true;
