
#pragma once

#include <sstream>
#include <string>
#include "primitives.hpp"

template <int precision=8>
auto serialize(const point_t &p) noexcept {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << p.x << "," << p.y << "," << p.z;
    return stream.str();
}
