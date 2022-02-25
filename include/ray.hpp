
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <glm/glm.hpp>

template <typename T>
struct ray_t {
    glm::vec<3,T> p,d{};
};
