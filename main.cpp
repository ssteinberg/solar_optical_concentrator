////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// standard library
#include <cmath>
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

// floating-point precision
using float_type = double;
constexpr int NUMBER_OF_DIMENSIONS = 2;
typedef linalg::vec<float_type, NUMBER_OF_DIMENSIONS> float_vec;
constexpr float_type A_LITTLE_BIT = 0.00001;

// global constants
constexpr float_type PI = 3.14159265358979;
constexpr float_type PI_OVER_TWO = PI / 2;
constexpr float_type DEG_TO_RAD = PI / 180;
constexpr float_type RAD_TO_DEG = 180 / PI;
constexpr float_type EPSILON = 0.01;
constexpr float_type EPSILON_OVER_TWO = EPSILON / 2;

// switches
constexpr bool INFINITE_SYSTEM = false;
constexpr bool FINITE_SYSTEM = !INFINITE_SYSTEM;
constexpr bool FLAT_SOURCE = true;
constexpr bool CYLINDRICAL_SOURCE = !FLAT_SOURCE;
constexpr bool FLAT_TARGET = true;
constexpr bool CYLINDRICAL_TARGET = !FLAT_TARGET;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// uniform rv c code
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
    Ray(const float_vec& o, const float_vec& d) : o(o), d(d) {}
};

// path
struct Path {

    std::vector<float_vec> vertices;

    void addVertex(const float_vec& v) {
        vertices.push_back(v);
    }

    void writePath(std::ofstream& file) {
        const int pathLength = vertices.size() - 1;
        if (pathLength > 0) {
            for (int i = 0; i < pathLength; ++i) {
                const float_vec v1 = vertices[i];
                const float_vec v2 = vertices[i + 1];
                file << v1.x << "," << v1.y << "," << v2.x << "," << v2.y << "\n";
            }
        }
    }

    void writeFinalSegment(std::ofstream& file) {
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
    CYLINDRICAL
};

// type
enum struct Type {
    SOURCE,
    TARGET,
    MIRROR,
    MIRROR_1,
    MIRROR_2,
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
    virtual bool intersect(const Ray& ray, HitInfo& hitInfo) = 0;
    virtual Ray sampleMeanRay() = 0;
    virtual Ray sampleDiffuseRay() = 0;
};

// line segment
struct LineSegment : Geometry {

    Shape shape = Shape::FLAT;
    Type type;
    float_vec p1, p2, n1, n2;

    LineSegment(const Type& t, const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) : type(t), p1(p1), p2(p2), n1(n1), n2(n2) {}

    bool intersect(const Ray& ray, HitInfo& hitInfo) override {
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

    Ray sampleMeanRay() override {
        const float_type s = PCG32::rand();
        const float_vec p = (1 - s) * p1 + s * p2;
        const float_vec n = (1 - s) * n1 + s * n2;
        return Ray(p + A_LITTLE_BIT * n, n);
    }

    Ray sampleDiffuseRay() override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(meanRay.o, rotDir);
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = -90; degrees <= 90; ++degrees) {
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

    void writeToFile(std::ofstream& file) {
        file << p1.x << "," << p1.y << "," << p2.x << "," << p2.y << "\n";
    }

};

// circle
struct Circle : Geometry {

    Shape shape = Shape::CYLINDRICAL;
    Type type;
    float_vec c;
    float_type r;

    Circle(const Type& t, const float_vec& c, const float_type& r) : type(t), c(c), r(r) {}

    bool intersect(const Ray& ray, HitInfo& hitInfo) override {
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

    Ray sampleMeanRay() override {
        const float_type theta = PCG32::rand() * 2 * PI;
        const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
        const float_vec n = (p - c) / r;
        return Ray(p, n);
    }

    Ray sampleDiffuseRay() override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return Ray(meanRay.o, rotDir);
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeRays() {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = 0; degrees <= 360; degrees += 2) {
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

};

// mirror
struct Mirror : Geometry {

    std::vector<LineSegment> segments;

    void addSegment(const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) {
        segments.emplace_back(Type::MIRROR, p1, p2, n1, n2);
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo) override {
        for (auto s : segments) {
            if (s.intersect(ray, hitInfo)) return true;
        }
        return false;
    }

    Ray sampleMeanRay() override {
        return segments[PCG32::rand() * segments.size()].sampleMeanRay();
    }

    Ray sampleDiffuseRay() override {
        return segments[PCG32::rand() * segments.size()].sampleDiffuseRay();
    }

    void writeToFile(std::ofstream& file) {
        for (auto s : segments) {
            s.writeToFile(file);
        }
    }

};

// two-mirror concentrator
struct TwoMirrorConcentrator : Geometry {

    Mirror m1, m2, barrier;

    void buildInfinite(const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max) {

        // initial conditions
        float_vec p(-L, 0), pp(f, 0), n(1, 0), np(-1, 0);
        float_type R(L + f), r(f);

        // numerical integration
        float_type B(0);
        while (std::abs(B) < B_max) {

            // step in beta
            const float_type B_new = B + dB;

            // precomputation
            const float_type sinB = std::sin(B);
            const float_type cosB = std::cos(B);
            const float_type d = R - 2 * (L + f);

            // (7)
            const float_type I = 1;
            const float_type S = 1;
            const float_type dy = (inv ? -1 : 1) * S / I * dB;
            const float_type dx = dy * (p.y - r * sinB) / (r * (cosB - 1) + 2 * (L + f));
            const float_vec p_new(p.x + dx, p.y + dy);

            // (11)
            const float_type dr = r * ((r + d) * sinB - p.y * cosB) / ((r + d) * cosB + p.y * sinB - r - R) * dB;
            const float_type r_new = r + dr;
            const float_vec pp_new = r_new * float_vec(std::cos(B_new), std::sin(B_new));

            // new directions
            float_vec K_int = pp_new - p_new;
            const float_type R_new = length(K_int);
            // const float_type R_new = 2 * (L + f) + p_new.x - r_new;
            K_int /= R_new;
            const float_vec K_out = -pp_new / r_new;
            const float_vec n_new = normalize(K_int - K_in);
            const float_vec np_new = normalize(K_out - K_int);

            // store coordinates and normals
            m1.addSegment(p, p_new, n, n_new);
            m2.addSegment(pp, pp_new, np, np_new);

            // update
            p = p_new;
            pp = pp_new;
            n = n_new;
            np = np_new;
            R = R_new;
            r = r_new;
            B = B_new;
        }
    }

    void buildFinite(const auto& Sa, const auto& SB, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max) {

        // initial conditions
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
            const float_type dB = Sa(a) * da / SB(B);
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
            m1.addSegment(p1, p1_new, n1, n1_new);
            m2.addSegment(p2, p2_new, n2, n2_new);
            m1.addSegment(float_vec(p1.x, -p1.y), float_vec(p1_new.x, -p1_new.y), float_vec(n1.x, -n1.y), float_vec(n1_new.x, -n1_new.y));
            m2.addSegment(float_vec(p2.x, -p2.y), float_vec(p2_new.x, -p2_new.y), float_vec(n2.x, -n2.y), float_vec(n2_new.x, -n2_new.y));

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
    }

    void buildBarrier() {
        const float_type x_max = 2 * std::max(std::abs(m1.segments.front().p1.x), std::abs(m2.segments.front().p1.x));
        const float_type y_max = 2 * std::max(std::abs(m1.segments.back().p2.y), std::abs(m2.segments.back().p2.y));
        const float_vec bottomRight(x_max, -y_max), topRight(x_max, y_max), topLeft(-x_max, y_max), bottomLeft(-x_max, -y_max);
        barrier.addSegment(bottomRight, topRight, float_vec(-1, 0), float_vec(-1, 0));
        barrier.addSegment(topRight, topLeft, float_vec(0, -1), float_vec(0, -1));
        barrier.addSegment(topRight, bottomRight, float_vec(1, 0), float_vec(1, 0));
        barrier.addSegment(bottomLeft, bottomRight, float_vec(0, 1), float_vec(0, 1));
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo) override {
        if (m1.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_1;
            return true;
        }
        if (m2.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_2;
            return true;
        }
        if (barrier.intersect(ray, hitInfo)) {
            hitInfo.t = Type::BARRIER;
            return true;
        }
        return false;
    }

    Ray sampleMeanRay() override {
        if (PCG32::rand() < 0.5) return m1.sampleMeanRay();
        else return m2.sampleMeanRay();
    }

    Ray sampleDiffuseRay() override {
        if (PCG32::rand() < 0.5) return m1.sampleDiffuseRay();
        else return m2.sampleDiffuseRay();
    }

};

// system
struct Design {

    std::vector<Geometry*> geometries;

    void addGeometry(Geometry* const geometry) {
        geometries.push_back(geometry);
    }

    bool intersect(const Ray& ray, HitInfo& minHitInfo) {
        bool hit = false;
        HitInfo tempMinHitInfo;
        minHitInfo.l = std::numeric_limits<float_type>::max();
        for (auto geometry : geometries) {
            if ((*geometry).intersect(ray, tempMinHitInfo)) {
                if (tempMinHitInfo.l < minHitInfo.l) {
                    hit = true;
                    minHitInfo = tempMinHitInfo;
                }
            }
        }
        return hit;
    }

    Path traceRay(const Ray& ray) {
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        int i = 0;
        while (i < 3) {
            if (HitInfo h; intersect(r, h)) {

                if (i == 0 && h.t != Type::MIRROR_1) return Path();
                if (i == 1 && h.t != Type::MIRROR_2) return Path();

                path.addVertex(h.p);
                float_vec reflectedDirection = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + A_LITTLE_BIT * reflectedDirection, reflectedDirection);
            }
            ++i;
        }
        return path;
    }

    std::vector<Path> rayTrace(const std::vector<Ray>& rays) {
        std::vector<Path> paths;
        for (auto ray : rays) {
            paths.push_back(traceRay(ray));
        }
        return paths;
    }

};



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



// main
int main(int argc, char* argv[]) {

    // finite system
    if (FINITE_SYSTEM) {

        // input
        const float_type f1(1), L(4), f2(1), da(0.001), a_max(110 * DEG_TO_RAD);
        const float_type radius(f1 / 100);

        // source
        LineSegment flatSource(Type::SOURCE, float_vec(-L, -radius), float_vec(-L, radius), float_vec(-1, 0), float_vec(-1, 0));
        Circle cylindricalSource(Type::SOURCE, float_vec(-L, 0), radius);
        const auto& Sa = [=](const float_type& alpha) {
            if (FLAT_SOURCE) return std::cos(alpha);
            else if (CYLINDRICAL_SOURCE) return float_type(1);
            else return float_type(0);
        };

        // target
        LineSegment flatTarget(Type::TARGET, float_vec(0, -radius), float_vec(0, radius), float_vec(1, 0), float_vec(1, 0));
        Circle cylindricalTarget(Type::TARGET, float_vec(0, 0), radius);
        const auto& SB = [=](const float_type& beta) {
            if (FLAT_TARGET) return std::cos(beta);
            else if (CYLINDRICAL_TARGET) return float_type(1);
            else return float_type(0);
        };

        // two-mirror concentrator
        TwoMirrorConcentrator tmc;
        tmc.buildFinite(Sa, SB, f1, L, f2, da, a_max);
        
        // design
        Design design;
        // if (FLAT_SOURCE) design.addGeometry(&flatSource);
        // else if (CYLINDRICAL_SOURCE) design.addGeometry(&cylindricalSource);
        // if (FLAT_TARGET) design.addGeometry(&flatTarget);
        // else if (CYLINDRICAL_TARGET) design.addGeometry(&cylindricalTarget);
        design.addGeometry(&tmc);

        // ray trace
        std::vector<Ray> rays;
        if (FLAT_SOURCE) {
            for (int i = 0; i < 100; ++i) {
                rays.push_back(flatSource.sampleMeanRay());
                // rays.push_back(flatSource.sampleDiffuseRay());
            }
        }
        if (CYLINDRICAL_SOURCE) {
            for (int i = 0; i < 100; ++i) {
                rays.push_back(cylindricalSource.sampleMeanRay());
                // rays.push_back(cylindricalSource.sampleDiffuseRay());
            }
        }
        std::vector<Path> paths = design.rayTrace(rays);

        // output
        std::ofstream file;
        file.open("data_2mc_fin_cyl_input.csv");
        file << FLAT_SOURCE << "," << FLAT_TARGET << "," << radius << "," << f1 << "," << L << "," << f2 << "," << da << "," << a_max << "\n";
        file.close();
        file.open("data_2mc_fin_cyl_m1.csv");
        tmc.m1.writeToFile(file);
        file.close();
        file.open("data_2mc_fin_cyl_m2.csv");
        tmc.m2.writeToFile(file);
        file.close();
        file.open("data_2mc_fin_cyl_paths.csv");
        for (auto path : paths) {
            path.writePath(file);
        }
        file.close();



        // extreme rays
        tmc.buildBarrier();
        if (FLAT_SOURCE) {
            auto extrRays = flatSource.generateExtremeRays();
            auto p1Paths = design.rayTrace(extrRays.first);
            auto p2Paths = design.rayTrace(extrRays.second);
            file.open("data_2mc_fin_extr_1.csv");
            for (auto path : p1Paths) path.writeFinalSegment(file);
            file.close();
            file.open("data_2mc_fin_extr_2.csv");
            for (auto path : p2Paths) path.writeFinalSegment(file);
            file.close();
        }
        if (CYLINDRICAL_SOURCE) {
            auto extrRays = cylindricalSource.generateExtremeRays();
            auto p1Paths = design.rayTrace(extrRays.first);
            auto p2Paths = design.rayTrace(extrRays.second);
            file.open("data_2mc_fin_extr_1.csv");
            for (auto path : p1Paths) path.writeFinalSegment(file);
            file.close();
            file.open("data_2mc_fin_extr_2.csv");
            for (auto path : p2Paths) path.writeFinalSegment(file);
            file.close();
        }

    }

}



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/*

// scale-invariant
if (SCALE_INVARIANT) {
    const TargetShape ts(TargetShape::CYLINDRICAL);
    const bool inv(true);
    constexpr float_type scale(10);
    constexpr float_type f(scale * 0.5);
    constexpr float_type L(scale * 3);
    constexpr float_type B_max(150 * DEG_TO_RAD);
    buildInf2MC(ts, inv, L, f, K_in, dB, B_max, "data_scale_inv.csv");
}

*/

/*

// parabolic mirror
// std::vector<std::tuple<float_vec, float_vec, float_vec>> PMa, PMb, PMc;
// float_type h = -1;
// while (h < 1) {
//     const Ray meanRay(float_vec(0, h), float_vec(-1, 0));
//     PM.traceRay(PMa, meanRay);
//     const Ray posExtrRay(float_vec(0, h), posExtrRayDir);
//     PM.traceRay(PMb, posExtrRay);
//     const Ray negExtrRay(float_vec(0, h), negExtrRayDir);
//     PM.traceRay(PMc, negExtrRay);
//     h += 0.01;
// }
// file.open("data_pm_a.csv");
// for (auto path : PMa) {
//     file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
//     file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
//     file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
// }
// file.close();
// file.open("data_pm_b.csv");
// for (auto path : PMb) {
//     file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
//     file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
//     file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
// }
// file.close();
// file.open("data_pm_c.csv");
// for (auto path : PMc) {
//     file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
//     file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
//     file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
// }
// file.close();

*/

/*

// // check intersection with mirror 1
// if (HitInfo hit1; M1.intersect(hit1, ray1)) {

//     // check intersection with mirror 2
//     const Ray ray2(hit1.P, normalize(ray1.d - 2 * dot(hit1.N, ray1.d) * hit1.N));
//     if (HitInfo hit2; M2.intersect(hit2, ray2)) {

//         // compute point through target and store path of ray
//         const float_type len = 6;
//         rayPaths.emplace_back(ray1.o, hit1.P, hit2.P, hit2.P + len * normalize(ray2.d - 2 * dot(hit2.N, ray2.d) * hit2.N));
//     }
// }



// // trace multiple rays across y direction
// void traceRaysIncrementally(const float_vec& rayDirection, const float_type& h_min, const float_type& h_inc, const float_type& h_max, const std::string& fileName) {

//     // trace
//     std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> paths;
//     float_type h = h_min;
//     while (h < h_max) {
//         const Ray ray(float_vec(0, h), rayDirection);
//         traceRay(paths, ray);
//         h += h_inc;
//     }

//     // write
//     std::ofstream file;
//     file.open(fileName);
//     for (auto rayPath : paths) {
//         file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
//         file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
//         file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
//         file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
//     }
//     file.close();
// }

*/

/*

// test sampling
std::ofstream file;
file.open("data_2mc_fin_cyl_paths.csv");
for (int i = 0; i < 100; ++i) {
    auto sample = source.sampleSurface();
    file << sample.x << "," << sample.y << "\n";
}
file.close();

# test sampling
plt.scatter(data_paths[:, 0], data_paths[:, 1])

*/

/*

// check all line segments in mirror
int segmentsSize = static_cast<int>(segments.size());
int normalsSize = static_cast<int>(normals.size());
int n = segmentsSize == normalsSize ? segmentsSize : 0;
for (int i = 0; i < n; ++i) {

    // check parallel
    float_vec segmentDirection = segments[i].second - segments[i].first;
    float_type denom = cross(segmentDirection, ray.d);
    if (denom == 0) continue;

    // intersection
    float_type t = cross((ray.o - segments[i].first), ray.d) / denom;
    if (0 - A_LITTLE_BIT <= t && t <= 1 + A_LITTLE_BIT) {
        hitInfo.p = segments[i].first + t * segmentDirection;
        hitInfo.n = normalize(t * normals[i].second + (1 - t) * normals[i].first);
        hitInfo.l = length(hitInfo.p - ray.o);
        return true;
    }
}

// no intersection
return false;

*/

/*

// std::cout << "ray.o: (" << r.o.x << ", " << r.o.y << ")" << std::endl;
// std::cout << "ray.d: (" << r.d.x << ", " << r.d.y << ")" << std::endl;
// std::cout << "hit.l: " << h.l << std::endl;
// std::cout << "hit.p: (" << h.p.x << ", " << h.p.y << ")" << std::endl;
// std::cout << "hit.n: (" << h.n.x << ", " << h.n.y << ")" << std::endl;
// std::cout << std::endl;

*/

/*

// check parallel
const float_vec segDir = seg.second - seg.first;
const float_type rayDir_x_segDir = cross(ray.d, segDir);
if (rayDir_x_segDir == 0) continue;

// intersection
const float_type s = cross((seg.first - ray.o), ray.d) / rayDir_x_segDir;
const float_type t = cross((seg.first - ray.o), segDir) / rayDir_x_segDir;
if (0 <= s && s <= 1 && 0 < t) {
    hitInfo.l = t;
    hitInfo.p = ray.o + hitInfo.l * ray.d;
    hitInfo.n = normalize((1 - s) * norms.first + s * norms.second);
    return true;
}

*/

/*

// infinite system
if (INFINITE_SYSTEM) {

    // input
    const float_vec K_in(-1, 0);
    const float_type dB(0.01);

    // flat target
    TwoMirrorConcentrator infFlat2MC;
    ParabolicMirror infFlatPM;
    if (FLAT_TARGET) {
        const TargetShape ts(TargetShape::FLAT);
        const bool inv(false);
        const float_type f(0.2);
        const float_type L(8 * f);
        const float_type B_max(80 * DEG_TO_RAD);
        infFlat2MC = buildInf2MC(ts, inv, L, f, K_in, dB, B_max, "data_2mc_inf_flat.csv");

        // parabolic mirror
        if (PARABOLIC_MIRROR) {
            const float_type dy(dB);
            infFlatPM = buildPM(f, dy, B_max, "data_pm_inf_flat.csv");
        }
    }

    // cylindrical target
    TwoMirrorConcentrator infCyl2MC;
    ParabolicMirror infCylPM;
    if (CYLINDRICAL_TARGET) {
        const TargetShape ts(TargetShape::CYLINDRICAL);
        const bool inv(true);
        const float_type f(0.5);
        const float_type L(6 * f);
        const float_type B_max(150 * DEG_TO_RAD);
        infCyl2MC = buildInf2MC(ts, inv, L, f, K_in, dB, B_max, "data_2mc_inf_cyl.csv");

        // parabolic mirror
        if (PARABOLIC_MIRROR) {
            const float_type dy(dB);
            infCylPM = buildPM(f, dy, B_max, "data_pm_inf_cyl.csv");
        }
    }

    // ray trace
    const float_type h_inc = 0.01;
    if (FLAT_TARGET) {
        infFlat2MC.traceRaysIncrementally(K_in, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_mean.csv");
        float_vec extrPosDir(-std::cos(EPSILON_OVER_TWO), std::sin(EPSILON_OVER_TWO));
        infFlat2MC.traceRaysIncrementally(extrPosDir, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_extr_pos.csv");
        float_vec extrNegDir(-std::cos(EPSILON_OVER_TWO), -std::sin(EPSILON_OVER_TWO));
        infFlat2MC.traceRaysIncrementally(extrNegDir, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_extr_neg.csv");
    }
    if (CYLINDRICAL_TARGET) {
        infCyl2MC.traceRaysIncrementally(K_in, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_mean.csv");
        float_vec extrPosDir(-std::cos(EPSILON_OVER_TWO), std::sin(EPSILON_OVER_TWO));
        infCyl2MC.traceRaysIncrementally(extrPosDir, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_extr_pos.csv");
        float_vec extrNegDir(-std::cos(EPSILON_OVER_TWO), -std::sin(EPSILON_OVER_TWO));
        infCyl2MC.traceRaysIncrementally(extrNegDir, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_extr_neg.csv");
    }
}

*/

/*

// build infinite two-mirror concentrator
TwoMirrorConcentrator buildInf2MC(const Shape& s, const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max,
    const std::string& fileName) {

    // build
    TwoMirrorConcentrator inf2MC;
    inf2MC.buildInfinite(s, inv, L, f, K_in, dB, B_max);

    // write
    std::ofstream file;
    file.open(fileName);

    // for (int i = 0, n = static_cast<int>(inf2MC.M1.segments.size()); i < n; ++i) {
    //     file << inf2MC.M1.segments[i].first.x << "," << inf2MC.M1.segments[i].first.y << ",";
    //     file << inf2MC.M2.segments[i].first.x << "," << inf2MC.M2.segments[i].first.y << "\n";
    // }
    // file << inf2MC.M1.segments.back().second.x << "," << inf2MC.M1.segments.back().second.y << ",";
    // file << inf2MC.M2.segments.back().second.x << "," << inf2MC.M2.segments.back().second.y << "\n";

    file.close();

    // return
    return inf2MC;
}

// build parabolic mirror
ParabolicMirror buildPM(const float_type& f, const float_type& dy, const float_type& B_max,
    const std::string& fileName) {

    // build
    ParabolicMirror PM;
    PM.build(f, dy, B_max);

    // write
    std::ofstream file;
    file.open(fileName);

    // for (auto segment : PM.M.segments) file << segment.first.x << "," << segment.first.y << "\n";

    file.close();

    // return
    return PM;
}

*/

/*

// parabolic mirror
struct ParabolicMirror : Geometry {

    Mirror M;

    void build(const float_type& f, const float_type& dy, const float_type& B_max) {
        float_vec p(-f, 0), n(1, 0);
        float_type B(0);
        while (B < B_max) {
            const float_type y_new = p.y + dy;
            const float_type x_new = y_new * y_new / (4 * f) - f;
            const float_vec p_new(x_new, y_new);
            const float_vec n_new = normalize(float_vec(2 * f, -p_new.y));
            M.addSegment(p, p_new, n, n_new);
            p = p_new;
            n = n_new;
            B = std::atan(p_new.y / p_new.x) * RAD_TO_DEG;
        }
    }

    void traceRay(std::vector<std::tuple<float_vec, float_vec, float_vec>>& rayPaths, const Ray& ray) {
        if (HitInfo hitInfo; M.intersect(ray, hitInfo)) {
            rayPaths.emplace_back(ray.o, hitInfo.p, hitInfo.p + normalize(ray.d - 2 * dot(hitInfo.n, ray.d) * hitInfo.n));
        }
    }

};

*/
