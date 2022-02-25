
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <limits>
#include <optional>
#include <glm/glm.hpp>

#include "aabb.hpp"
#include "ray.hpp"

#include "primitives.hpp"


template <typename T>
struct intersection_t {
    T d{ std::numeric_limits<T>::infinity() };
    glm::vec<3,T> n,refl{};
};

// aabb-aabb test
template <typename T>
bool intersect(const aabb<T>& o1, const aabb<T>& o2) noexcept {
    return glm::all((glm::lessThanEqual(o1.a, o2.a) && 
                     glm::lessThan(o2.a, o1.b)) ||
                    (glm::lessThanEqual(o2.a, o1.a) && 
                     glm::lessThan(o1.a, o2.b)));
}

// ray-aabb test
template <typename T>
bool intersect(const ray_t<T> &r, const aabb<T> &bb) noexcept {
    auto inv_d = glm::vec<3,T>{1,1,1} / r.d;
    const auto &t1 = (bb.a - r.p) * inv_d;
    const auto &t2 = (bb.b - r.p) * inv_d;

    const auto &tmin3 = glm::min(t1, t2);
    const auto &tmax3 = glm::max(t1, t2);
    const auto &tmin  = glm::max(tmin3.x, glm::max(tmin3.y, tmin3.z));
    const auto &tmax  = glm::min(tmax3.x, glm::min(tmax3.y, tmax3.z));

    return tmax>=tmin && tmax>.0f;
}

inline std::optional<intersection_t<f_t>> intersect_ray_line(const ray_t<f_t> &r, const line_t &l, bool back=false) {
    const auto &v1=l.a, &v2=l.b;
    const auto v = v2-v1;
    if (v==point_t{0,0,0})
        return std::nullopt;

    const auto D = -r.d.x*(v1.y-v2.y) + r.d.y*(v1.x-v2.x);
    if (glm::abs(D)<1e-20)
        return std::nullopt;
    const auto t = ((r.p.x-v1.x)*(v1.y-v2.y) - (r.p.y-v1.y)*(v1.x-v2.x)) / D;
    const auto u = (r.d.x*(r.p.y-v1.y) - r.d.y*(r.p.x-v1.x)) / D;
    if (t<=1e-8 || u<0 || u>1)
        return std::nullopt;

    intersection_t<f_t> intr;
    intr.d=t;

    const auto& n = glm::mix(l.n1,l.n2,u);
    if (!back && glm::dot(n,r.d)>=0)
        return std::nullopt;
    intr.n = glm::normalize(n);
    intr.refl = glm::reflect(r.d,intr.n);
    
    return intr;
}

inline std::optional<intersection_t<f_t>> intersect_ray_tri(const ray_t<f_t> &r, const tri_t &l, bool back=false) {
    assert(false && "todo");
    return std::nullopt;
}

template <typename T>
bool solve_quadratic(T a, T b, T c, T &x0, T &x1) { 
    auto d = b * b - 4 * a * c; 
    if (d < 0) return false; 
    else if (d == 0) x0 = x1 = -0.5 * b/a;
    else { 
        auto q = -0.5 * (b + (b>0?T(1):T(-1))*glm::sqrt(d));
        x0 = q/a;
        x1 = c/q;
        if (x0 > x1) std::swap(x0, x1);
    }
 
    return true; 
} 
inline std::optional<intersection_t<f_t>> intersect_ray_sphere(const ray_t<f_t> &ray, const ball_t &ball, bool front=true) {
    const auto& L = ray.p - ball.c;
    const auto  a = glm::dot(ray.d,ray.d); 
    const auto  b = 2*glm::dot(ray.d,L); 
    const auto  c = glm::dot(L,L) - sqr(ball.r);
    f_t t0,t1;
    if (!solve_quadratic(a, b, c, t0, t1)) return std::nullopt;
    if (t0<0 && t1<0) return std::nullopt;
    if (t0<0) t0=t1;

    intersection_t<f_t> intr;
    intr.d = front?t0:t1;

    const auto p = ray.p+intr.d*ray.d;
    assert(glm::abs(glm::length(p-ball.c)-ball.r)<1e-12);
    intr.n = (p-ball.c)/ball.r;
    intr.refl = glm::reflect(ray.d,intr.n);

    return intr;
}

inline auto closest_ray_point(const ray_t<f_t> &ray, const point_t &p) {
    const auto t = glm::max<f_t>(0, glm::dot(ray.d,p-ray.p));
    return ray.p + ray.d * t;
}
