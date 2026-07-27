////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// standard library
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <ranges>
#include <tuple>
#include <vector>

// linear algebra
#include "linalg.h"
using namespace linalg::aliases;

// parallelization
#include <omp.h>

// floating-point precision
using float_type = double;
typedef linalg::vec<float_type, 2> float_vec;

// global constants
constexpr float_type PI = 3.14159265358979;
constexpr float_type PI_OVER_TWO = PI / 2;
constexpr float_type ONE_OVER_PI = 1 / PI;
constexpr float_type DEG_TO_RAD = PI / 180;
constexpr float_type RAD_TO_DEG = 180 / PI;
constexpr float_type EPSILON = 0.05;
constexpr float_type EPSILON_OVER_TWO = EPSILON / 2;

// source
constexpr bool FLAT_SOURCE = false;
constexpr bool CYL_SOURCE = false;
constexpr bool ELL_SOURCE = true;

// target
constexpr bool FLAT_TARGET = false;
constexpr bool CYL_TARGET = false;
constexpr bool ELL_TARGET = true;

// switches
constexpr bool IGNORE_SOURCE = true;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// uniform random variable
#include <stdint.h>
namespace PCG32 {
	static uint64_t mcg_state = 0xcafef00dd15ea5e5u; // must be odd
	static uint64_t const multiplier = 6364136223846793005u;
	uint32_t pcg32_fast(void) {
		uint64_t x = mcg_state;
		const unsigned count = (unsigned)(x >> 61);
		mcg_state = x * multiplier;
		x ^= x >> 22;
		return (uint32_t)(x >> (22 + count));
	}
	float rand() {
		return float(double(pcg32_fast()) / 4294967296.0);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// ray
struct Ray {
    float_vec o, d;
    Ray() : o(), d(float_vec(0.0, 1.0)) {}
    Ray(const float_vec& o, const float_vec& d) : o(o), d(d) {}
};

// path
struct Path {

    std::vector<float_vec> vertices;

    void addVertex(const float_vec& v) {
        vertices.push_back(v);
    }

    void writePath(std::ofstream& file) const {
        const int pathLength = vertices.size() - 1;
        if (pathLength > 0) {
            for (int i = 0; i < pathLength; ++i) {
                const float_vec v1 = vertices[i];
                const float_vec v2 = vertices[i + 1];
                file << v1.x << "," << v1.y << "," << v2.x << "," << v2.y << "\n";
            }
        }
    }

    void writeFinalSegment(std::ofstream& file) const {
        const int pathLength = vertices.size() - 1;
        if (pathLength > 0) {
            const float_vec v1 = vertices[pathLength - 1];
            const float_vec v2 = vertices[pathLength];
            file << v1.x << "," << v1.y << "," << v2.x << "," << v2.y << "\n";
        }
    }
    
};

// shape
enum struct Shape {
    FLAT,
    CYLINDRICAL,
    ELLIPTICAL,
    PARABOLIC,
    CONSTRUCTED
};

// type
enum struct Type {
    SOURCE,
    TARGET,
    CONCENTRATOR,
    MIRROR1,
    MIRROR2,
    BARRIER
};

// hit info
struct HitInfo {
    float_type l;
    float_vec p, n;
    Type t;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// geometry
struct Geometry {
    Shape shape;
    Type type;
    Geometry(const Shape& s, const Type& t) : shape(s), type(t) {}
    virtual float_vec getCentre() const = 0;
    virtual float_type getLength() const = 0;
    virtual float_type getEllipseSemiA() const = 0;
    virtual float_type getEllipseSemiB() const = 0;
    virtual bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const = 0;
    virtual Ray sampleSourceRay() const = 0;
    virtual Ray sampleMeanRay() const = 0;
    virtual std::pair<Ray, float_type> sampleDiffuseRay() const = 0;
    virtual std::vector<Ray> generatePointRays() const = 0;
    virtual std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const = 0;
    virtual std::vector<Ray> generateFinalPlotRays() const = 0;
};

// line segment
struct LineSegment : Geometry {

    float_vec p1, p2, n1, n2;
    LineSegment(const Type& t, const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) : Geometry(Shape::FLAT, t), p1(p1), p2(p2), n1(n1), n2(n2) {}
    float_vec getCentre() const override { return (p1 + p2) / 2; }
    float_type getLength() const override { return length(p1 - p2); }
    float_type getEllipseSemiA() const override { return 0; }
    float_type getEllipseSemiB() const override { return 0; }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        for (const auto& t : types) {
            if (t == type) {
                if (IGNORE_SOURCE && type == Type::SOURCE) return false;
                else if (dot(n1, -ray.d) <= 0 && type == Type::TARGET ) return false;
                else {
                    const float_vec segDir = p2 - p1;
                    const float_type rayDir_x_segDir = cross(ray.d, segDir);
                    if (rayDir_x_segDir == 0) return false;
                    const float_type s = cross((p1 - ray.o), ray.d) / rayDir_x_segDir;
                    const float_type t = cross((p1 - ray.o), segDir) / rayDir_x_segDir;
                    if (0 <= s && s <= 1 && 0 < t) {
                        hitInfo.l = t;
                        hitInfo.p = ray.o + hitInfo.l * ray.d;
                        hitInfo.n = normalize((1 - s) * n1 + s * n2);
                        hitInfo.t = type;
                        return true;
                    }
                    return false;
                }
            }
        }
        return false;
    }

    Ray sampleSourceRay() const override {
        const float_vec n(-1, 0);
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(getCentre() + float_vec(0, 0.001), rotDir); // can change ray.o to check imaging
    }

    Ray sampleMeanRay() const override {
        const float_type s = PCG32::rand();
        const float_vec p = (1 - s) * p1 + s * p2;
        const float_vec n = (1 - s) * n1 + s * n2;
        return Ray(p, n);
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return std::pair(Ray(meanRay.o, rotDir), theta);
    }

    std::vector<Ray> generatePointRays() const override {
        const Ray r = sampleMeanRay();
        std::vector<Ray> pointRays;
        for (int degrees = -90; degrees <= 90; degrees += 5) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec rotDir(r.d.x * cos0 - r.d.y * sin0, r.d.x * sin0 + r.d.y * cos0);
            pointRays.emplace_back(r.o, rotDir);
        }
        return pointRays;
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = -90; degrees <= 90; degrees += 5) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec p1RotDir(n1.x * cos0 - n1.y * sin0, n1.x * sin0 + n1.y * cos0);
            const float_vec p2RotDir(n2.x * cos0 - n2.y * sin0, n2.x * sin0 + n2.y * cos0);
            p1ExtrRays.emplace_back(p1, p1RotDir);
            p2ExtrRays.emplace_back(p2, p2RotDir);
        }
        return std::pair(p1ExtrRays, p2ExtrRays);
    }

    std::vector<Ray> generateFinalPlotRays() const override {
        std::vector<Ray> rays;
        for (int degrees = -90; degrees <= 90; degrees += 15) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec p1RotDir(n1.x * cos0 - n1.y * sin0, n1.x * sin0 + n1.y * cos0);
            rays.emplace_back(getCentre(), p1RotDir);
        }
        return rays;
    }

    void writeLineSegment(std::ofstream& file) const {
        file << p1.x << "," << p1.y << "," << p2.x << "," << p2.y << "\n";
    }

};

// circle
struct Circle : Geometry {

    float_vec c; float_type r;
    Circle(const Type& t, const float_vec& c, const float_type& r) : Geometry(Shape::CYLINDRICAL, t), c(c), r(r) {}
    float_vec getCentre() const override { return c; }
    float_type getLength() const override { return PI * 2 * r; }
    float_type getEllipseSemiA() const override { return 0; }
    float_type getEllipseSemiB() const override { return 0; }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        for (const auto& t : types) {
            if (t == type) {
                if (IGNORE_SOURCE && type == Type::SOURCE) return false;
                else {
                    const float_vec oc = c - ray.o;
                    const float_type hyp2 = dot(oc, oc);
                    if (hyp2 < r * r) return false;
                    const float_type hyp = std::sqrt(hyp2);
                    const float_vec ocDir = oc / hyp;
                    const float_type cos0 = dot(ray.d, ocDir);
                    if (cos0 <= 0) return false;
                    const float_type adj = hyp * cos0;
                    const float_type d2 = hyp2 - adj * adj;
                    if (d2 > r * r) return false;
                    if (d2 == r * r) {
                        hitInfo.l = adj;
                        hitInfo.p = ray.o + hitInfo.l * ray.d;
                    }
                    else if (d2 < r * r) {
                        hitInfo.l = adj - std::sqrt(r * r - d2);
                        hitInfo.p = ray.o + hitInfo.l * ray.d;
                    }
                    hitInfo.n = normalize(hitInfo.p - c);
                    hitInfo.t = type;
                    return true;
                }
            }
        }
        return false;
    }

    Ray sampleSourceRay() const override {
        const float_vec n(-1, 0);
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(getCentre(), rotDir);
    }

    Ray sampleMeanRay() const override {
        const float_type theta = PCG32::rand() * 2 * PI;
        const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
        const float_vec n = (p - c) / r;
        return Ray(p, n);
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return std::pair(Ray(meanRay.o, rotDir), theta * RAD_TO_DEG);
    }

    std::vector<Ray> generatePointRays() const override {
        const Ray r = sampleMeanRay();
        std::vector<Ray> pointRays;
        for (int degrees = -90; degrees <= 90; degrees += 5) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec rotDir(r.d.x * cos0 - r.d.y * sin0, r.d.x * sin0 + r.d.y * cos0);
            pointRays.emplace_back(r.o, rotDir);
        }
        return pointRays;
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = 0; degrees < 360; degrees += 5) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
            const float_vec n = (p - c) / r;
            const float_vec p1RotDir(-n.y, n.x);
            const float_vec p2RotDir(n.y, -n.x);
            p1ExtrRays.emplace_back(p, p1RotDir);
            p2ExtrRays.emplace_back(p, p2RotDir);
        }
        return std::pair(p1ExtrRays, p2ExtrRays);
    }

    std::vector<Ray> generateFinalPlotRays() const override {
        std::vector<Ray> rays;
        for (int degrees = -170; degrees <= 170; degrees += 20) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec n = float_vec(-1, 0);
            const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
            rays.emplace_back(getCentre() * rotDir, rotDir);
        }
        return rays;
    }

};

// parabola
struct Parabola : Geometry {

    float_type a, c;
    Parabola(const Type& t, const float_type& a, const float_type& c) : Geometry(Shape::PARABOLIC, t), a(a), c(c) {}
    float_vec getCentre() const override { return float_vec(c, 0); }
    float_type getLength() const override { return float_type(0); }
    float_type getEllipseSemiA() const override { return 0; }
    float_type getEllipseSemiB() const override { return 0; }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        for (const auto& t : types) {
            if (t == type) {
                if (IGNORE_SOURCE && type == Type::SOURCE) return false;
                else {
                    if (ray.d.x == -1 || ray.d.x == -1 || ray.d.y == 0) {
                        const float_type s = (a * ray.o.y * ray.o.y - ray.o.x + c) / ray.d.x;
                        if (s < 0) return false;
                        hitInfo.l = s;
                        hitInfo.p = ray.o + hitInfo.l * ray.d;
                        hitInfo.n = normalize(float_vec(-1, 2 * a * hitInfo.p.y));
                        hitInfo.t = type;
                        return true;
                    }
                    const float_type A = a * ray.d.y * ray.d.y;
                    const float_type B = 2 * a * ray.o.y * ray.d.y - ray.d.x;
                    const float_type C = a * ray.o.y * ray.o.y + c - ray.o.x;
                    const float_type D = B * B - 4 * A * C;
                    if (D <= 0 || A == 0) return false;
                    const float_type invA = 1 / A;
                    const float_type sqrtD = std::sqrt(D);
                    const float_type s1 = 0.5 * (-B - std::copysign(1, B) * sqrtD) * invA;
                    const float_type s2 = C * invA / s1;
                    float_type s;
                    if (s1 >= 0 && s2 >= 0) s = std::min(s1, s2);
                    else if (s1 >= 0) s = s1;
                    else if (s2 >= 0) s = s2;
                    else return false;
                    hitInfo.l = s;
                    hitInfo.p = ray.o + hitInfo.l * ray.d;
                    if (hitInfo.p.y > 0) hitInfo.n = normalize(float_vec(1, -2 * a * hitInfo.p.y));
                    else hitInfo.n = normalize(float_vec(-1, 2 * a * hitInfo.p.y));
                    hitInfo.t = type;
                    return true;
                }
                return false;
            }
        }
        return false;
    }

    Ray sampleSourceRay() const override {
        const float_vec n(-1, 0);
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(getCentre(), rotDir);
    }

    Ray sampleMeanRay() const override {
        return Ray(float_vec(0, 0), float_vec(0, 0));
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        return std::pair(Ray(float_vec(0, 0), float_vec(0, 0)), float_type(0));
    }

    std::vector<Ray> generatePointRays() const override {
        return std::vector<Ray>();
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    std::vector<Ray> generateFinalPlotRays() const override {
        return std::vector<Ray>();
    }

};

// ellipse
struct Ellipse : Geometry {

    float_type a, b;
    float_vec c;
    Ellipse(const Type& t, const float_type& a, const float_type& b, const float_vec& c) : Geometry(Shape::ELLIPTICAL, t), a(a), b(b), c(c) {}
    float_vec getCentre() const override { return c; }
    float_type getLength() const override { return std::max(a, b); }
    float_type getEllipseSemiA() const override { return a; }
    float_type getEllipseSemiB() const override { return b; }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        for (const auto& t : types) {
            if (t == type) {
                if (IGNORE_SOURCE && type == Type::SOURCE) return false;
                else {
                    const float_vec oc = ray.o - c;
                    const float_type inva2 = 1 / (a * a);
                    const float_type invb2 = 1 / (b * b);
                    const float_type A = ray.d.x * ray.d.x * inva2 + ray.d.y * ray.d.y * invb2;
                    const float_type B = 2 * (oc.x * ray.d.x * inva2 + oc.y * ray.d.y * invb2);
                    const float_type C = oc.x * oc.x * inva2 + oc.y * oc.y * invb2 - 1;
                    const float_type D = B * B - 4 * A * C;
                    if (D < 0) return false;
                    const float_type invA = 1 / A;
                    const float_type sqrtD = std::sqrt(D);
                    const float_type s1 = 0.5 * (-B - std::copysign(1, B) * sqrtD) * invA;
                    const float_type s2 = C * invA / s1;
                    float_type s;
                    if (s1 >= 0 && s2 >= 0) s = std::min(s1, s2);
                    else if (s1 >= 0 || s2 >= 0) s = std::max(s1, s2); // added to fix interior hit detection
                    else return false;
                    hitInfo.l = s;
                    hitInfo.p = ray.o + hitInfo.l * ray.d;
                    hitInfo.n = normalize(float_vec((hitInfo.p.x - c.x) * inva2, (hitInfo.p.y - c.y) * invb2));
                    hitInfo.t = type;
                    return true;
                }
            }
        }
        return false;
    }

    Ray sampleSourceRay() const override { return Ray(getCentre(), float_vec(0, 1)); }

    Ray sampleMeanRay() const override {
        float_type theta(PCG32::rand() * 2 * PI), costheta(std::cos(theta)), sintheta(std::sin(theta));
        while (PCG32::rand() > std::sqrt(b * b * costheta * costheta + a * a * sintheta * sintheta) / std::max(a, b)) {
            theta = PCG32::rand() * 2 * PI;
            costheta = std::cos(theta);
            sintheta = std::sin(theta);
        }
        const float_vec p = c + float_vec(a * costheta, b * sintheta);
        const float_vec n = normalize(float_vec(b * costheta, a * sintheta));
        return Ray(p, n);
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type costheta = std::cos(theta);
        const float_type sintheta = std::sin(theta);
        const float_vec rotDir(n.x * costheta - n.y * sintheta, n.x * sintheta + n.y * costheta);
        return std::pair(Ray(meanRay.o, rotDir), theta);
    }

    std::vector<Ray> generatePointRays() const override { return std::vector<Ray>(); }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        std::vector<Ray> extrRays1, extrRays2;
        for (int degrees = 0; degrees < 360; degrees += 5) {
            const float_type theta(degrees * DEG_TO_RAD), costheta(std::cos(theta)), sintheta(std::sin(theta));
            const float_vec p = c + float_vec(a * costheta, b * sintheta);
            const float_vec n = normalize(float_vec(b * costheta, a * sintheta));
            const float_vec left(-n.y, n.x);
            const float_vec right(n.y, -n.x);
            extrRays1.emplace_back(p, left);
            extrRays2.emplace_back(p, right);
        }
        return std::pair(extrRays1, extrRays2);
    }

    std::vector<Ray> generateFinalPlotRays() const override { return std::vector<Ray>(); }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// mirror
struct Mirror : Geometry {

    std::vector<LineSegment> segments;
    Mirror(const Type& t) : Geometry(Shape::CONSTRUCTED, t) {}
    void reset() { segments.clear(); }
    float_vec getCentre() const override { return float_vec(0, 0); }

    float_type getLength() const override {
        float_type l(0);
        for (const auto& s : segments) l += s.getLength();
        return l;
    }

    float_type getEllipseSemiA() const override { return 0; }
    float_type getEllipseSemiB() const override { return 0; }

    void addSegment(const Type& t, const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) {
        segments.emplace_back(t, p1, p2, n1, n2);
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        for (const auto& t : types) {
            const std::vector<Type> temp{t};
            if (t == type) {
                for (auto s : segments) if (s.intersect(ray, hitInfo, temp)) return true;
                return false;
            }
        }
        return false;
    }

    Ray sampleSourceRay() const override {
        const float_vec n(-1, 0);
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(getCentre(), rotDir);
    }

    Ray sampleMeanRay() const override {
        return segments[PCG32::rand() * segments.size()].sampleMeanRay();
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        return segments[PCG32::rand() * segments.size()].sampleDiffuseRay();
    }

    std::vector<Ray> generatePointRays() const override {
        return std::vector<Ray>();
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    std::vector<Ray> generateFinalPlotRays() const override {
        return std::vector<Ray>();
    }

    void writeMirror(std::ofstream& file) const {
        for (auto s : segments) {
            s.writeLineSegment(file);
        }
    }

};

// two-mirror concentrator
struct TwoMirrorConcentrator : Geometry {

    Mirror m1a = Mirror(Type::MIRROR1);
    Mirror m1b = Mirror(Type::MIRROR1);
    Mirror m2a = Mirror(Type::MIRROR2);
    Mirror m2b = Mirror(Type::MIRROR2);
    Mirror barrier = Mirror(Type::BARRIER);

    TwoMirrorConcentrator() : Geometry(Shape::CONSTRUCTED, Type::CONCENTRATOR) {}

    void reset() {
        m1a.reset();
        m1b.reset();
        m2a.reset();
        m2b.reset();
        barrier.reset();
    }

    float_vec getCentre() const override {
        return float_vec(0, 0);
    }

    float_type getLength() const override {
        return m1a.getLength() + m1b.getLength() + m2a.getLength() + m2b.getLength();
    }

    float_type getEllipseSemiA() const override { return 0; }
    float_type getEllipseSemiB() const override { return 0; }

    float_type buildFin(const bool& inv, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max, const float_type& w, const float_type& hl) {

        // reset
        reset();

        // initial conditions
        const auto& Sa = [=](const float_type& alpha) {
            if (!ELL_SOURCE) {
                if (FLAT_SOURCE) return std::cos(alpha);
                else return float_type(1);
            }
            const float_type a(hl), b(hl * 2);
            const float_type cosalpha(std::cos(alpha)), sinalpha(std::sin(alpha));
            return 2 * std::sqrt(a * a * sinalpha * sinalpha + b * b * cosalpha * cosalpha);
        };
        const auto& SB = [=](const float_type& beta) {
            if (!ELL_TARGET) {
                if (FLAT_SOURCE && FLAT_TARGET) return std::cos(beta) * w;
                if (FLAT_SOURCE && !FLAT_TARGET) return ONE_OVER_PI;
                if (!FLAT_SOURCE && FLAT_TARGET) return PI * std::cos(beta);
                return float_type(1);
            }
            const float_type a(hl * 2), b(hl);
            const float_type cosbeta(std::cos(beta)), sinbeta(std::sin(beta));
            return 2 * std::sqrt(a * a * sinbeta * sinbeta + b * b * cosbeta * cosbeta);
        };
        float_type a(0), B(0), r1(f1), r2(f2);
        float_vec p1(-f1 - L, 0), p2(f2, 0), n1(1, 0), n2(-1, 0);
        const float_type F(2 * f1 + L + 2 * f2);
        const float_vec source(-L, 0), target(0, 0);

        // numerical integration
        while (a < a_max) {

            // precomputation
            const float_type cosa(std::cos(a)), sina(std::sin(a));
            const float_type cosB(std::cos(B)), sinB(std::sin(B));
            const float_type cosaB(std::cos(a + B)), sinaB(std::sin(a + B));

            // M1
            const float_type a_new = a + da;
            const float_type dr1 = (r1 * r2 * sinaB + L * r1 * sina) / (r2 * cosaB + L * cosa - r2 + F) * da;
            const float_type r1_new = r1 + dr1;
            const float_vec p1_new(-r1_new * std::cos(a_new) - L, r1_new * std::sin(a_new));

            // M2
            const float_type dB = (inv ? -1 : 1) * Sa(a) / SB(B) * da;
            const float_type B_new = B + dB;
            const float_type dr2 = (r1 * r2 * sinaB + L * r2 * sinB) / (r1 * cosaB + L * cosB - r1 + F) * dB;
            const float_type r2_new = r2 + dr2;
            const float_vec p2_new(r2_new * std::cos(B_new), r2_new * std::sin(B_new));

            // ray directions
            const float_vec v1((p1_new - source) / r1_new);
            const float_vec u((p2_new - p1_new) / (F - r1_new - r2_new));
            const float_vec v2((target - p2_new) / r2_new);

            // normals
            const float_vec n1_new(normalize(u - v1));
            const float_vec n2_new(normalize(v2 - u));

            // store coordinates and normals
            m1a.addSegment(m1a.type, p1, p1_new, n1, n1_new);
            m1b.addSegment(m1b.type, float_vec(p1.x, -p1.y), float_vec(p1_new.x, -p1_new.y), float_vec(n1.x, -n1.y), float_vec(n1_new.x, -n1_new.y));
            m2a.addSegment(m2a.type, p2, p2_new, n2, n2_new);
            m2b.addSegment(m2b.type, float_vec(p2.x, -p2.y), float_vec(p2_new.x, -p2_new.y), float_vec(n2.x, -n2.y), float_vec(n2_new.x, -n2_new.y));

            // update
            a = a_new;
            r1 = r1_new;
            p1 = p1_new;
            B = B_new;
            r2 = r2_new;
            p2 = p2_new;
            n1 = n1_new;
            n2 = n2_new;
        }

        // build barrier
        buildBarrier();

        // final beta
        return B;
    }

    void buildBarrier() {
        const float_type x_max = 1.25 * std::max(std::abs(m1a.segments.front().p1.x), std::abs(m2a.segments.front().p1.x));
        const float_type y_max = 1.25 * std::max(std::abs(m1a.segments.back().p2.y), std::abs(m2a.segments.back().p2.y));
        const float_vec bottomRight(x_max, -y_max), topRight(x_max, y_max), topLeft(-x_max, y_max), bottomLeft(-x_max, -y_max);
        barrier.addSegment(barrier.type, bottomRight, topRight, float_vec(-1, 0), float_vec(-1, 0));
        barrier.addSegment(barrier.type, topRight, topLeft, float_vec(0, -1), float_vec(0, -1));
        barrier.addSegment(barrier.type, topLeft, bottomLeft, float_vec(1, 0), float_vec(1, 0));
        barrier.addSegment(barrier.type, bottomLeft, bottomRight, float_vec(0, 1), float_vec(0, 1));
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo, const std::vector<Type>& types) const override {
        // bool hit = false;
        // HitInfo tempHitInfo;
        // hitInfo.l = std::numeric_limits<float_type>::max();
        // for (const auto& t : types) {
        //     const std::vector<Type> temp{t};
        //     if (t == Type::MIRROR1) {
        //         if (m1a.intersect(ray, tempHitInfo, temp) && tempHitInfo.l < hitInfo.l) {
        //             hit = true;
        //             hitInfo = tempHitInfo;
        //         }
        //         if (m1b.intersect(ray, tempHitInfo, temp) && tempHitInfo.l < hitInfo.l) {
        //             hit = true;
        //             hitInfo = tempHitInfo;
        //         }
        //     }
        //     if (t == Type::MIRROR2) {
        //         if (m2a.intersect(ray, tempHitInfo, temp) && tempHitInfo.l < hitInfo.l) {
        //             hit = true;
        //             hitInfo = tempHitInfo;
        //         }
        //         if (m2b.intersect(ray, tempHitInfo, temp) && tempHitInfo.l < hitInfo.l) {
        //             hit = true;
        //             hitInfo = tempHitInfo;
        //         }
        //     }
        //     if (t == Type::BARRIER) {
        //         if (barrier.intersect(ray, tempHitInfo, temp) && tempHitInfo.l < hitInfo.l) {
        //             hit = true;
        //             hitInfo = tempHitInfo;
        //         }
        //     }
        // }
        // return hit;

        for (const auto& t : types) {
            const std::vector<Type> temp{t};
            if (t == Type::MIRROR1) {
                if (m1a.intersect(ray, hitInfo, temp)) return true;
                if (m1b.intersect(ray, hitInfo, temp)) return true;
            }
            if (t == Type::MIRROR2) {
                if (m2a.intersect(ray, hitInfo, temp)) return true;
                if (m2b.intersect(ray, hitInfo, temp)) return true;
            }
            if (t == Type::BARRIER) if (barrier.intersect(ray, hitInfo, temp)) return true;
        }
        return false;
    }

    Ray sampleSourceRay() const override {
        const float_vec n(-1, 0);
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(getCentre(), rotDir);
    }

    Ray sampleMeanRay() const override {
        const float_type Randy = PCG32::rand();
        if (Randy < 0.25) return m1a.sampleMeanRay();
        if (0.25 <= Randy && Randy < 0.5) return m1b.sampleMeanRay();
        if (0.5 <= Randy && Randy < 0.75) return m2a.sampleMeanRay();
        return m2b.sampleMeanRay();
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const float_type Randolf = PCG32::rand();
        if (Randolf < 0.25) return m1a.sampleDiffuseRay();
        if (0.25 <= Randolf && Randolf < 0.5) return m1b.sampleDiffuseRay();
        if (0.5 <= Randolf && Randolf < 0.75) return m2a.sampleDiffuseRay();
        return m2b.sampleDiffuseRay();
    }

    std::vector<Ray> generatePointRays() const override {
        return std::vector<Ray>();
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    std::vector<Ray> generateFinalPlotRays() const override {
        return std::vector<Ray>();
    }

    void writeTwoMirrorConcentrator(const std::string& filePath) const {
        std::ofstream file;
        file.open(filePath + "mirror1a.csv");
        m1a.writeMirror(file);
        file.close();
        file.open(filePath + "mirror1b.csv");
        m1b.writeMirror(file);
        file.close();
        file.open(filePath + "mirror2a.csv");
        m2a.writeMirror(file);
        file.close();
        file.open(filePath + "mirror2b.csv");
        m2b.writeMirror(file);
        file.close();
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// design
struct Design {

    Geometry* source;
    Geometry* target;
    std::vector<Geometry*> geometries;

    void addGeometry(Geometry* const g) {
        geometries.push_back(g);
        if (g->type == Type::SOURCE) source = g;
        else if (g->type == Type::TARGET) target = g;
    }

    bool intersect(const Ray& ray, HitInfo& minHitInfo, const std::vector<Type>& types) const {
        bool hit = false;
        HitInfo tempMinHitInfo;
        minHitInfo.l = std::numeric_limits<float_type>::max();
        for (const auto& geometry : geometries) {
            if ((*geometry).intersect(ray, tempMinHitInfo, types)) {
                if (tempMinHitInfo.l < minHitInfo.l) {
                    hit = true;
                    minHitInfo = tempMinHitInfo;
                }
            }
        }
        return hit;
    }

    void traceSourceRays(const std::string& filePath, const int& numRays) {
        std::vector<std::vector<Path>> threadVectors(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            Path path;
            const Ray r1 = source->sampleSourceRay();
            path.addVertex(r1.o);
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    path.addVertex(h2.p);
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::TARGET};
                    if (HitInfo h3; intersect(r3, h3, temp3)) path.addVertex(h3.p);
                }
            }
            threadVectors[threadID].push_back(path);
        }
        std::ofstream file;
        std::vector<Path> paths;
        for (const auto& v : threadVectors) paths.insert(paths.end(), v.begin(), v.end());
        file.open(filePath + "source.csv");
        for (const auto& p : paths) p.writePath(file);
        file.close();
    }

    void traceDiffuseRays(const std::string& filePath, const int& numRays) {
        std::vector<std::vector<Path>> threadVectors(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            Path path;
            const Ray r1 = source->sampleDiffuseRay().first;
            path.addVertex(r1.o);
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    path.addVertex(h2.p);
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::TARGET};
                    if (HitInfo h3; intersect(r3, h3, temp3)) path.addVertex(h3.p);
                }
            }
            threadVectors[threadID].push_back(path);
        }
        std::ofstream file;
        std::vector<Path> paths;
        for (const auto& v : threadVectors) paths.insert(paths.end(), v.begin(), v.end());
        file.open(filePath + "diffuse.csv");
        for (const auto& p : paths) p.writePath(file);
        file.close();
    }

    void tracePointRays(const std::string& filePath, const int& numPoints) const {
        for (int i = 0; i < numPoints; ++i) {
            auto pointRays = source->generatePointRays();
            std::vector<std::vector<Path>> threadVectors(omp_get_max_threads());
            #pragma omp parallel for schedule(dynamic, 1)
            for (const auto& ray : pointRays) {
                int threadID = omp_get_thread_num();
                Path path;
                path.addVertex(ray.o);
                const Ray r1 = ray;
                const std::vector<Type> temp1{Type::MIRROR1};
                if (HitInfo h1; intersect(r1, h1, temp1)) {
                    path.addVertex(h1.p);
                    const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                    const std::vector<Type> temp2{Type::MIRROR2};
                    if (HitInfo h2; intersect(r2, h2, temp2)) {
                        path.addVertex(h2.p);
                        const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                        const std::vector<Type> temp3{Type::BARRIER};
                        if (HitInfo h3; intersect(r3, h3, temp3)) path.addVertex(h3.p);
                    }
                }
                threadVectors[threadID].push_back(path);
            }
            std::vector<Path> paths;
            for (const auto& v : threadVectors) paths.insert(paths.end(), v.begin(), v.end());
            std::ofstream file;
            file.open(filePath + "pointrays" + std::to_string(i) + ".csv");
            for (const auto& path : paths) path.writePath(file);
            file.close();
        }
    }

    void traceExtremeRays(const std::string& filePath) const {
        auto extremeRays = source->generateExtremeRays();
        std::vector<std::vector<Path>> threadVectors1(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (const auto& ray : extremeRays.first) {
            int threadID = omp_get_thread_num();
            Path path;
            path.addVertex(ray.o);
            const Ray r1 = ray;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    path.addVertex(h2.p);
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::BARRIER};
                    if (HitInfo h3; intersect(r3, h3, temp3)) path.addVertex(h3.p);
                }
            }
            threadVectors1[threadID].push_back(path);
        }
        std::vector<Path> paths1;
        for (const auto& v : threadVectors1) paths1.insert(paths1.end(), v.begin(), v.end());
        std::vector<std::vector<Path>> threadVectors2(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (const auto& ray : extremeRays.second) {
            int threadID = omp_get_thread_num();
            Path path;
            path.addVertex(ray.o);
            const Ray r1 = ray;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    path.addVertex(h2.p);
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::BARRIER};
                    if (HitInfo h3; intersect(r3, h3, temp3)) path.addVertex(h3.p);
                }
            }
            threadVectors2[threadID].push_back(path);
        }
        std::vector<Path> paths2;
        for (const auto& v : threadVectors2) paths2.insert(paths2.end(), v.begin(), v.end());
        std::ofstream file;
        file.open(filePath + "extreme1.csv");
        for (const auto& path : paths1) path.writePath(file);
        file.close();
        file.open(filePath + "extreme2.csv");
        for (const auto& path : paths2) path.writePath(file);
        file.close();
    }

    void traceExtremeEllipse(const std::string& filePath) const {
        auto extremeRays = source->generateExtremeRays();
        std::vector<std::vector<Path>> threadVectors1(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (const auto& ray : extremeRays.first) {
            int threadID = omp_get_thread_num();
            Path path;
            path.addVertex(ray.o);
            const Ray r1 = ray;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::BARRIER};
                if (HitInfo h2; intersect(r2, h2, temp2)) path.addVertex(h2.p);
            }
            threadVectors1[threadID].push_back(path);
        }
        std::vector<Path> paths1;
        for (const auto& v : threadVectors1) paths1.insert(paths1.end(), v.begin(), v.end());
        std::vector<std::vector<Path>> threadVectors2(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (const auto& ray : extremeRays.second) {
            int threadID = omp_get_thread_num();
            Path path;
            path.addVertex(ray.o);
            const Ray r1 = ray;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                path.addVertex(h1.p);
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::BARRIER};
                if (HitInfo h2; intersect(r2, h2, temp2)) path.addVertex(h2.p);
            }
            threadVectors2[threadID].push_back(path);
        }
        std::vector<Path> paths2;
        for (const auto& v : threadVectors2) paths2.insert(paths2.end(), v.begin(), v.end());
        std::ofstream file;
        file.open(filePath + "extreme1.csv");
        for (const auto& path : paths1) path.writeFinalSegment(file);
        file.close();
        file.open(filePath + "extreme2.csv");
        for (const auto& path : paths2) path.writeFinalSegment(file);
        file.close();
    }

    void tracePhaseSpace(const std::string& filePath, const int& numRays) const {
        std::vector<std::vector<std::vector<float_type>>> threadVectors(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            auto sample = source->sampleDiffuseRay();
            const Ray r1 = sample.first;
            const float_type source_angle = sample.second;
            const float_vec source_pos = r1.o - source->getCentre();
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::TARGET};
                    if (HitInfo h3; intersect(r3, h3, temp3)) {
                        if (source->shape == Shape::ELLIPTICAL && target->shape == Shape::ELLIPTICAL) {
                            float_type source_theta = std::atan2(source_pos.y * source->getEllipseSemiA(), -source_pos.x * source->getEllipseSemiB());
                            if (source_theta < 0) source_theta += 2 * PI;
                            float_type target_theta = std::atan2(h3.p.y * target->getEllipseSemiA(), h3.p.x * target->getEllipseSemiB());
                            if (target_theta < 0) target_theta += 2 * PI;
                            std::vector<float_type> temp = {cross(-r3.d, h3.n), target_theta * RAD_TO_DEG, std::sin(source_angle), source_theta * RAD_TO_DEG};
                            threadVectors[threadID].emplace_back(std::move(temp));
                        }
                    }
                }
            }
        }
        std::ofstream file;
        std::vector<std::vector<float_type>> data;
        for (const auto& v : threadVectors) data.insert(data.end(), v.begin(), v.end());
        file.open(filePath + "phase.csv");
        for (const auto& d : data) file << d[0] << "," << d[1] << "," << d[2] << "," << d[3] << "\n";
        file.close();
    }

    void tracePhaseEllipse(const std::string& filePath, const int& numRays) const {
        std::vector<std::vector<std::vector<float_type>>> threadVectors(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            auto sample = source->sampleDiffuseRay();
            const Ray r1 = sample.first;
            const float_type source_angle = sample.second;
            const float_vec source_pos = r1.o - source->getCentre();
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::TARGET};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    if (source->shape == Shape::FLAT && target->shape == Shape::FLAT) {
                        std::vector<float_type> temp = {cross(-r2.d, h2.n), h2.p.y, std::sin(source_angle), source_pos.y};
                        threadVectors[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::CYLINDRICAL && target->shape == Shape::CYLINDRICAL) {
                        float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                        if (source_theta < 0) source_theta += 2 * PI;
                        float_type target_theta = std::atan2(h2.p.y, h2.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r2.d, h2.n), target_theta * RAD_TO_DEG, std::sin(source_angle), source_theta * RAD_TO_DEG};
                        threadVectors[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::FLAT && target->shape == Shape::CYLINDRICAL) {
                        float_type target_theta = std::atan2(h2.p.y, h2.p.x);
                        if (target_theta < 0) target_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r2.d, h2.n), target_theta * RAD_TO_DEG, std::sin(source_angle), source_pos.y};
                        threadVectors[threadID].emplace_back(std::move(temp));
                    }
                    else if (source->shape == Shape::CYLINDRICAL && target->shape == Shape::FLAT) {
                        float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                        if (source_theta < 0) source_theta += 2 * PI;
                        std::vector<float_type> temp = {cross(-r2.d, h2.n), h2.p.y, std::sin(source_angle), source_theta * RAD_TO_DEG};
                        threadVectors[threadID].emplace_back(std::move(temp));
                    }
                }
            }
        }
        std::ofstream file;
        std::vector<std::vector<float_type>> data;
        for (const auto& v : threadVectors) data.insert(data.end(), v.begin(), v.end());
        file.open(filePath + "phase.csv");
        for (const auto& d : data) file << d[0] << "," << d[1] << "," << d[2] << "," << d[3] << "\n";
        file.close();
    }

    float_type traceHitData(const int& numRays) const {
        std::vector<int> threadCounts(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            const Ray r1 = source->sampleDiffuseRay().first;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::MIRROR2};
                if (HitInfo h2; intersect(r2, h2, temp2)) {
                    const Ray r3 = Ray(h2.p, normalize(r2.d - 2 * dot(r2.d, h2.n) * h2.n));
                    const std::vector<Type> temp3{Type::TARGET};
                    if (HitInfo h3; intersect(r3, h3, temp3)) threadCounts[threadID] += 1;
                }
            }

        }
        int finalCount = 0;
        for (const auto& count : threadCounts) finalCount += count;
        return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
    }

    float_type traceHitEllipse(const int& numRays) const {
        std::vector<int> threadCounts(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            const Ray r1 = source->sampleDiffuseRay().first;
            const std::vector<Type> temp1{Type::MIRROR1};
            if (HitInfo h1; intersect(r1, h1, temp1)) {
                const Ray r2 = Ray(h1.p, normalize(r1.d - 2 * dot(r1.d, h1.n) * h1.n));
                const std::vector<Type> temp2{Type::TARGET};
                if (HitInfo h2; intersect(r2, h2, temp2)) threadCounts[threadID] += 1;
            }
        }
        int finalCount = 0;
        for (const auto& count : threadCounts) finalCount += count;
        return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



// main
int main(int argc, char* argv[]) {



    // TWO-MIRROR CONCENTRATOR
    if (true) {

        // output
        std::string outputDataPath = "data2mc/";
        if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
        std::ofstream file;

        // input
        const bool inv(false);
        const float_type f1(1), L(30), f2(1), da(0.000001), a_max(158.45 * DEG_TO_RAD), w(1), hl(f1 / 1000);

        // two-mirror concentrator
        TwoMirrorConcentrator tmc;
        const float_type B_max = tmc.buildFin(inv, f1, L, f2, da, a_max, w, hl);
        tmc.writeTwoMirrorConcentrator(outputDataPath);

        // write parameters
        file.open(outputDataPath + "input.csv");
        file << inv << "," << hl << "," << f1 << "," << L << "," << f2 << "," << da << "," << a_max * RAD_TO_DEG << "," << w << "," << B_max * RAD_TO_DEG << "\n";
        file.close();

        // design
        Design design;
        design.addGeometry(&tmc);
        Ellipse ellSrc(Type::SOURCE, hl, hl * 2, float_vec(-L, 0));
        Ellipse ellTgt(Type::TARGET, hl * 2, hl, float_vec(0, 0));
        design.addGeometry(&ellSrc);
        design.addGeometry(&ellTgt);

        // hits
        const int numberOfDesigns = 20;
        std::vector<Design> designs(numberOfDesigns);
        std::vector<Ellipse> ellTgts;
        float_type increment = 2 * hl / numberOfDesigns;
        for (int i = 0; i < numberOfDesigns; ++i) {
            ellTgts.emplace_back(Type::TARGET, (i + 1) * increment * 2, (i + 1) * increment, float_vec(0, 0));
        }
        for (int i = 0; i < numberOfDesigns; ++i) {
            designs[i].addGeometry(&ellSrc);
            designs[i].addGeometry(&ellTgts[i]);
            designs[i].addGeometry(&tmc);
        }
        std::vector<float_vec> hitData;
        const int numberOfTrials = 1;
        const int numberOfRays = 10000;
        for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitData(numberOfRays));
        file.open(outputDataPath + "hits.csv");
        for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
        file.close();
    }



    // ELLIPTICAL CONCENTRATOR
    if (false) {

        // output
        std::string outputDataPath = "dataell/";
        if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
        std::ofstream file;

        // input
        const float_type f1(1), L(30), hl(f1 / 1000);
        file.open(outputDataPath + "input.csv");
        file << FLAT_SOURCE << "," << FLAT_TARGET << "," << hl << "," << f1 << "," << L << "\n";
        file.close();

        // elliptical concentrator
        Ellipse ec(Type::MIRROR1, f1 + L / 2, std::sqrt(f1 * (f1 + L)), float_vec(-L / 2, 0));
        Ellipse ellSrc(Type::SOURCE, hl, hl * 2, float_vec(-L, 0));
        Ellipse ellTgt(Type::TARGET, hl * 2, hl, float_vec(0, 0));

        // barrier
        const float_type size = (2 * f1 + L) * 2;
        const float_vec bottomRight(size, -size), topRight(size, size), topLeft(-size, size), bottomLeft(-size, -size);
        LineSegment rightBarrier(Type::BARRIER, bottomRight, topRight, float_vec(-1, 0), float_vec(-1, 0));
        LineSegment topBarrier(Type::BARRIER, topRight, topLeft, float_vec(0, -1), float_vec(0, -1));
        LineSegment leftBarrier(Type::BARRIER, topLeft, bottomLeft, float_vec(1, 0), float_vec(1, 0));
        LineSegment bottomBarrier(Type::BARRIER, bottomLeft, bottomRight, float_vec(0, 1), float_vec(0, 1));

        // design
        Design design;
        design.addGeometry(&ec);
        design.addGeometry(&ellSrc);
        design.addGeometry(&ellTgt);
        design.addGeometry(&rightBarrier);
        design.addGeometry(&topBarrier);
        design.addGeometry(&leftBarrier);
        design.addGeometry(&bottomBarrier);

        // hits
        const int numberOfDesigns = 20;
        std::vector<Design> designs(numberOfDesigns);
        std::vector<Ellipse> ellTgts;
        float_type increment = 2 * hl / numberOfDesigns;
        for (int i = 0; i < numberOfDesigns; ++i) {
            ellTgts.emplace_back(Type::TARGET, (i + 1) * increment * 2, (i + 1) * increment, float_vec(0, 0));
        }
        for (int i = 0; i < numberOfDesigns; ++i) {
            designs[i].addGeometry(&ellSrc);
            designs[i].addGeometry(&ellTgts[i]);
            designs[i].addGeometry(&ec);
        }
        std::vector<float_vec> hitData;
        const int numberOfTrials = 1;
        const int numberOfRays = 100000;
        for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitEllipse(numberOfRays));
        file.open(outputDataPath + "hits.csv");
        for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
        file.close();
    }



}



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////