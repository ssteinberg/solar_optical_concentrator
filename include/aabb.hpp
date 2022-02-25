
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <glm/glm.hpp>


template <typename T>
struct aabb {
    glm::vec<3,T> a{},b{};
};

template <typename T> 
auto splice(const aabb<T>& a, const glm::bvec3 &b) noexcept {
    const auto& mid = (a.a+a.b)/T(2);
    return aabb<T>{ glm::mix(a.a,mid,b), glm::mix(mid,a.b,b) };
}
// Containment test
template <typename T>
bool operator<=(const aabb<T>& o1, const aabb<T>& o2) noexcept {
    return glm::all(glm::lessThanEqual(o1.b, o2.b)) && 
           glm::all(glm::lessThanEqual(o2.a, o1.a));
}
// aabb of aabbs
template <typename T>
auto operator+(const aabb<T>& o1, const aabb<T>& o2) noexcept {
    return aabb<T>{ glm::min(o1.a,o2.a), glm::max(o1.b,o2.b) };
}
template <typename T>
auto operator+=(aabb<T>& o1, const aabb<T>& o2) noexcept {
    o1 = o1+o2;
    return o1;
}
