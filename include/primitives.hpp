
#pragma once

#include <glm/glm.hpp>
#include "common.hpp"
#include "aabb.hpp"


using f_t = double;

using point_t = glm::vec<3, f_t>;

enum class primitive_type : std::uint8_t {
    ball,
    line,
    tri,
};

template <typename T>
class primitive_base {
public:
    virtual aabb<T> AABB() const noexcept = 0;
};

struct primitive_t : primitive_base<f_t> {
    std::int8_t mirror{};
    bool is_target{};
    primitive_type type{};

    primitive_t(std::int8_t m, primitive_type type, bool is_target=false) noexcept : mirror(m),type(type),is_target(is_target) {}
    primitive_t(primitive_t &&) noexcept = default;
    primitive_t(const primitive_t &) noexcept = default;
    primitive_t& operator=(primitive_t &&) noexcept = default;
    primitive_t& operator=(const primitive_t &) noexcept = default;
};

struct line_t : primitive_t { 
    point_t a{}, b{};
    point_t n1{}, n2{};

    line_t(point_t a, point_t b, point_t n1, point_t n2, std::int8_t m, bool is_target=false) noexcept : primitive_t(m,primitive_type::line,is_target),a(a),b(b),n1(n1),n2(n2) {}
    line_t(line_t &&) noexcept = default;
    line_t(const line_t &) noexcept = default;
    line_t& operator=(line_t &&) noexcept = default;
    line_t& operator=(const line_t &) noexcept = default;
    
    aabb<f_t> AABB() const noexcept { return { glm::min(a,b), glm::max(a,b) }; }
};

struct tri_t : primitive_t { 
    point_t a{}, b{}, c{};
    point_t n1{};
    point_t n2{};
    point_t n3{};

    tri_t(point_t a, point_t b, point_t c, point_t n1, point_t n2, point_t n3, std::int8_t m, bool is_target=false) noexcept 
        : primitive_t(m,primitive_type::tri,is_target),a(a),b(b),c(c),n1(n1),n2(n2),n3(n3) {}
    tri_t(tri_t &&) noexcept = default;
    tri_t(const tri_t &) noexcept = default;
    tri_t& operator=(tri_t &&) noexcept = default;
    tri_t& operator=(const tri_t &) noexcept = default;

    aabb<f_t> AABB() const noexcept { return { glm::min(a,glm::min(b,c)), glm::max(a,glm::max(b,c)) }; }
};

struct ball_t : primitive_t {
    point_t c{};
    f_t r{};

    ball_t(point_t c, f_t r) noexcept : primitive_t(0,primitive_type::ball,true),c(c),r(r) {}
    ball_t(ball_t &&) noexcept = default;
    ball_t(const ball_t &) noexcept = default;
    ball_t& operator=(ball_t &&) noexcept = default;
    ball_t& operator=(const ball_t &) noexcept = default;

    aabb<f_t> AABB() const noexcept { return { c-point_t{r}, c+point_t{r} }; }
};

inline auto flip(const line_t& l) noexcept {
    line_t ret{ l };
    ret.a = l.b;
    ret.b = l.a;
    return ret;
}
inline auto flip(const tri_t& t) noexcept {
    tri_t ret{ t };
    ret.a = t.a;
    ret.b = t.c;
    ret.c = t.b;
    return ret;
}
