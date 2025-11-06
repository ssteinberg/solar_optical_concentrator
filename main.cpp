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

// doinking intersection points
constexpr bool DOINKING = true;
constexpr float_type DOINK = 1e-6;

// global constants
constexpr float_type PI = 3.14159265358979;
constexpr float_type PI_OVER_TWO = PI / 2;
constexpr float_type ONE_OVER_PI = 1 / PI;
constexpr float_type DEG_TO_RAD = PI / 180;
constexpr float_type RAD_TO_DEG = 180 / PI;
constexpr float_type EPSILON = 0.01;
constexpr float_type EPSILON_OVER_TWO = EPSILON / 2;

// finite system
constexpr bool BUILD_FINITE_SYSTEM = true;

// infinite system
constexpr bool BUILD_INFINITE_SYSTEM = false;
constexpr bool BUILD_INFINITE_ARBITRARY_SYSTEM = false;

// source and target
constexpr bool FLAT_SOURCE = true;
constexpr bool FLAT_TARGET = true;
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
    CONSTRUCTED
};

// type
enum struct Type {
    SOURCE,
    TARGET,
    MIRROR,
    MIRROR_1a,
    MIRROR_1b,
    MIRROR_2a,
    MIRROR_2b,
    BARRIER
};

// hit info
struct HitInfo {
    float_type l;
    float_vec p, n;
    Type t;
    float_type ml;
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
    virtual bool intersect(const Ray& ray, HitInfo& hitInfo) const = 0;
    virtual Ray sampleMeanRay() const = 0;
    virtual std::pair<Ray, float_type> sampleDiffuseRay() const = 0;
    virtual std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const = 0;
    virtual std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const = 0;
};

// line segment
struct LineSegment : Geometry {

    float_vec p1, p2, n1, n2;

    LineSegment(const Type& t, const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) : Geometry(Shape::FLAT, t), p1(p1), p2(p2), n1(n1), n2(n2) {}

    float_vec getCentre() const override { return (p1 + p2) / 2; }

    float_type getLength() const override { return length(p1 - p2); }

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override {
        if (IGNORE_SOURCE && type == Type::SOURCE) return false;
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

    Ray sampleMeanRay() const override {
        const float_type s = PCG32::rand();
        const float_vec p = (1 - s) * p1 + s * p2;
        const float_vec n = (1 - s) * n1 + s * n2;
        return Ray(p + (DOINKING ? DOINK * n : float_vec(0, 0)), n);
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return std::pair(Ray(meanRay.o, rotDir), std::abs(theta * RAD_TO_DEG));
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = -90; degrees <= 90; ++degrees) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_type cos0 = std::cos(theta);
            const float_type sin0 = std::sin(theta);
            const float_vec p1RotDir(n1.x * cos0 - n1.y * sin0, n1.x * sin0 + n1.y * cos0);
            const float_vec p2RotDir(n2.x * cos0 - n2.y * sin0, n2.x * sin0 + n2.y * cos0);
            p1ExtrRays.emplace_back(p1 + (DOINKING ? DOINK * n1 : float_vec(0, 0)), p1RotDir);
            p2ExtrRays.emplace_back(p2 + (DOINKING ? DOINK * n2 : float_vec(0, 0)), p2RotDir);
        }
        return std::pair(p1ExtrRays, p2ExtrRays);
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override {
        std::vector<Ray> posExtrRays, negExtrRays;
        const float_type inc(1 / static_cast<float_type>(numRays));
        float_type l(0);
        while (l <= 1) {
            const float_vec p((1 - l) * p1 + l * p2);
            const float_vec n((1 - l) * n1 + l * n2);
            const float_type cosPos(std::cos(EPSILON_OVER_TWO)), sinPos(std::sin(EPSILON_OVER_TWO));
            const float_type cosNeg(std::cos(-EPSILON_OVER_TWO)), sinNeg(std::sin(-EPSILON_OVER_TWO));
            const float_vec posRotDir(n.x * cosPos - n.y * sinPos, n.x * sinPos + n.y * cosPos);
            const float_vec negRotDir(n.x * cosNeg - n.y * sinNeg, n.x * sinNeg + n.y * cosNeg);
            posExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), posRotDir);
            negExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), negRotDir);
            l += inc;
        }
        return std::pair(posExtrRays, negExtrRays);
    }

    void writeLineSegment(std::ofstream& file) const {
        file << p1.x << "," << p1.y << "," << p2.x << "," << p2.y << "\n";
    }

};

// circle
struct Circle : Geometry {

    float_vec c;
    float_type r;

    Circle(const Type& t, const float_vec& c, const float_type& r) : Geometry(Shape::CYLINDRICAL, t), c(c), r(r) {}

    float_vec getCentre() const override { return c; }

    float_type getLength() const override { return PI * 2 * r; }

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override {
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

    Ray sampleMeanRay() const override {
        const float_type theta = PCG32::rand() * 2 * PI;
        const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
        const float_vec n = (p - c) / r;
        return Ray(p + (DOINKING ? DOINK * n : float_vec(0, 0)), n);
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        const Ray meanRay = sampleMeanRay();
        const float_vec n = meanRay.d;
        const float_type theta = std::asin(2 * PCG32::rand() - 1);
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        return std::pair(Ray(meanRay.o, rotDir), std::abs(theta * RAD_TO_DEG));
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override {
        std::vector<Ray> p1ExtrRays, p2ExtrRays;
        for (int degrees = 0; degrees <= 360; degrees += 2) {
            const float_type theta = degrees * DEG_TO_RAD;
            const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
            const float_vec n = (p - c) / r;
            const float_vec p1RotDir(-n.y, n.x);
            const float_vec p2RotDir(n.y, -n.x);
            p1ExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), p1RotDir);
            p2ExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), p2RotDir);
        }
        return std::pair(p1ExtrRays, p2ExtrRays);
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

};



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



// bounding box
struct BoundingBox {

	float_vec min, max;

    BoundingBox() : min(std::numeric_limits<float_type>::min(), std::numeric_limits<float_type>::min()),
        max(std::numeric_limits<float_type>::max(), std::numeric_limits<float_type>::max()) {}

    void fit(const LineSegment& l) {
        if (l.p1.x < min.x) min.x = l.p1.x;
        if (l.p1.y < min.y) min.y = l.p1.y;
        if (l.p2.x < min.x) min.x = l.p2.x;
        if (l.p2.y < min.y) min.y = l.p2.y;
        if (l.p1.x > max.x) max.x = l.p1.x;
        if (l.p1.y > max.y) max.y = l.p1.y;
        if (l.p2.x > max.x) max.x = l.p2.x;
        if (l.p2.y > max.y) max.y = l.p2.y;
    }

    bool intersect(const Ray& ray, HitInfo& minHit) const {
        float_type tx1 = (min.x - ray.o.x) / ray.d.x;
        float_type ty1 = (min.y - ray.o.y) / ray.d.y;
        float_type tx2 = (max.x - ray.o.x) / ray.d.x;
        float_type ty2 = (max.y - ray.o.y) / ray.d.y;
        if (tx1 > tx2) {
            const float_type temp = tx1;
            tx1 = tx2;
            tx2 = temp;
        }
        if (ty1 > ty2) {
            const float_type temp = ty1;
            ty1 = ty2;
            ty2 = temp;
        }
        float_type t1 = tx1;
        if (t1 < ty1) t1 = ty1;
        float_type t2 = tx2;
        if (t2 > ty2) t2 = ty2;
        if (t1 > t2) return false;
        if ((t1 < 0) && (t2 < 0)) return false;
        minHit.l = t1;
        return true;
    }

};

// node
struct Node {

    int size;
    BoundingBox bb;
    Node* left = nullptr;
    Node* right = nullptr;
    const LineSegment* ls = nullptr;

    Node(const std::vector<LineSegment>& v, const int& li, const int& ri) : size(ri - li + 1) {
        if (size > 1) {
            for (int i = li; i <= ri; ++i) bb.fit(v[i]);
            left = new Node(v, li, li + size / 2 - 1);
            right = new Node(v, li + size / 2, ri);
        }
        else ls = &v.at(li);
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo) const {
        if (size == 1) return ls->intersect(ray, hitInfo);
        if (bb.intersect(ray, hitInfo)) {
            if (left->intersect(ray, hitInfo)) return true;
            if (right->intersect(ray, hitInfo)) return true;
        }
        return false;
    }

};

// tree
struct Tree {
    Node* root = nullptr;
    Tree() = default;
    Tree(const std::vector<LineSegment>& v) { root = new Node(v, 0, static_cast<int>(v.size()) - 1); }
    bool intersect(const Ray& ray, HitInfo& hitInfo) const { if (root) return root->intersect(ray, hitInfo); return false; }
};



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



// mirror
struct Mirror : Geometry {

    std::vector<LineSegment> segments;
    Tree tree;

    Mirror() : Geometry(Shape::CONSTRUCTED, Type::MIRROR) {}

    void reset() { segments.clear(); }

    float_vec getCentre() const override { return float_vec(0, 0); }

    float_type getLength() const override {
        float_type l(0);
        for (const auto& s : segments) l += s.getLength();
        return l;
    }

    void addSegment(const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) {
        segments.emplace_back(Type::MIRROR, p1, p2, n1, n2);
    }

    void buildTree() {
        tree = Tree(segments);
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override {
        float_type length = 0;
        for (auto s : segments) {
            length += s.getLength();
            if (s.intersect(ray, hitInfo)) {
                hitInfo.ml = length;
                return true;
            }
        }
        return false;
        // return tree.intersect(ray, hitInfo);
    }

    Ray sampleMeanRay() const override {
        return segments[PCG32::rand() * segments.size()].sampleMeanRay();
    }

    std::pair<Ray, float_type> sampleDiffuseRay() const override {
        return segments[PCG32::rand() * segments.size()].sampleDiffuseRay();
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    void writeMirror(std::ofstream& file) const {
        for (auto s : segments) {
            s.writeLineSegment(file);
        }
    }

};



// two-mirror concentrator
struct TwoMirrorConcentrator : Geometry {

    Mirror m1a, m1b, m2a, m2b, barrier;

    TwoMirrorConcentrator() : Geometry(Shape::CONSTRUCTED, Type::MIRROR) {}

    void reset() {
        m1a.reset();
        m1b.reset();
        m2a.reset();
        m2b.reset();
        barrier.reset();
    }

    float_vec getCentre() const override { return float_vec(0, 0); }

    float_type getLength() const override {
        return m1a.getLength() + m1b.getLength() + m2a.getLength() + m2b.getLength();
    }

    void buildFin(const bool& inv, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max) {

        // reset
        reset();

        // initial conditions
        const auto& Sa = [=](const float_type& alpha) {
            if (FLAT_SOURCE) return std::cos(alpha);
            else return float_type(1);
        };
        const auto& SB = [=](const float_type& beta) {
            if (FLAT_SOURCE && FLAT_TARGET) return std::cos(beta);
            if (FLAT_SOURCE && !FLAT_TARGET) return ONE_OVER_PI;
            if (!FLAT_SOURCE && FLAT_TARGET) return PI * std::cos(beta);
            return float_type(1);
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
            m1a.addSegment(p1, p1_new, n1, n1_new);
            m1b.addSegment(float_vec(p1.x, -p1.y), float_vec(p1_new.x, -p1_new.y), float_vec(n1.x, -n1.y), float_vec(n1_new.x, -n1_new.y));
            m2a.addSegment(p2, p2_new, n2, n2_new);
            m2b.addSegment(float_vec(p2.x, -p2.y), float_vec(p2_new.x, -p2_new.y), float_vec(n2.x, -n2.y), float_vec(n2_new.x, -n2_new.y));

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

        // build trees
        // m1a.buildTree();
        // m1b.buildTree();
        // m2a.buildTree();
        // m2b.buildTree();

        // build barrier
        buildBarrier();
    }

    void buildInf(const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max) {

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
            const float_type I(1), S(FLAT_TARGET ? cosB : 1);
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
            m1a.addSegment(p, p_new, n, n_new);
            m1b.addSegment(float_vec(p.x, -p.y), float_vec(p_new.x, -p_new.y), float_vec(n.x, -n.y), float_vec(n_new.x, -n_new.y));
            m2a.addSegment(pp, pp_new, np, np_new);
            m2b.addSegment(float_vec(pp.x, -pp.y), float_vec(pp_new.x, -pp_new.y), float_vec(np.x, -np.y), float_vec(np_new.x, -np_new.y));

            // update
            p = p_new;
            pp = pp_new;
            n = n_new;
            np = np_new;
            R = R_new;
            r = r_new;
            B = B_new;
        }

        // build barrier
        buildBarrier();
    }



    std::pair<float_type, float_vec> calcCone(const float_vec& apex, const std::vector<float_vec>& vertices) {
        float_type sina_max(0);
        float_vec u(0, 0);
        for (const auto& v1 : vertices) {
            for (const auto& v2 : vertices) {
                const float_vec l1(normalize(v1 - apex)), l2(normalize(v2 - apex));
                const float_type sina(std::abs(cross(l1, l2)));
                if (sina > sina_max) {
                    sina_max = sina;
                    u = normalize(l1 + l2);
                }
            }
        }
        return std::pair(std::asin(sina_max), u);
    }

    bool traceRay(const Ray& r, HitInfo& h, Mirror& t) const {
        bool hit = false;
        HitInfo temp;
        h.l = std::numeric_limits<float_type>::max();
        if(t.intersect(r, temp)) {
            if (temp.l < h.l) {
                hit = true;
                h = temp;
            }
        }
        return hit;
    }

    void buildInfArb(const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max) {

        // reset
        reset();

        // flat target
        const float_vec tp1(0, -EPSILON_OVER_TWO), tp2(0, EPSILON_OVER_TWO);
        const float_vec tn(1, 0);
        Mirror t; t.addSegment(tp1, tp2, tn, tn);
        const std::vector<float_vec> tpts = {tp1, tp2};

        /*

        // triangular target
        // const float_vec tp1(1, 1), tp2(-2, 0), tp3(0, -1);
        // const float_vec tt1(tp1 - tp2), tt2(tp2 - tp3), tt3(tp3 - tp1);
        // const float_vec tn1(-tt1.y, tt1.x), tn2(-tt2.y, tt2.x), tn3(-tt3.y, tt3.x);
        // Mirror t;
        // t.addSegment(tp1, tp2, tn1, tn1);
        // t.addSegment(tp2, tp3, tn2, tn2);
        // t.addSegment(tp3, tp1, tn3, tn3);
        // const std::vector<float_vec> tpts = {tp1, tp2, tp3};

        */

        // initial conditions
        float_vec p1(-L, 0), p2(f, 0), n1(1, 0);
        float_type R(L + f), l(f);

        // calculate cone to initialize K_out and n2
        auto cone(calcCone(p2, tpts));
        float_vec K_out(cone.second);
        float_vec n2(normalize(p2 - p1) - K_out);

        // trace ray to initialize r
        if (HitInfo h; traceRay(Ray(p2, K_out), h, t)) {
            float_type a(cone.first), r(h.l);
            const float_type a_0(cone.first), F(2 * L + f + r);

            // numerical integration
            float_type B(0);
            while (std::abs(B) < B_max) {

                // step in beta
                const float_type B_new(B + dB);

                // precomputation
                const float_type sinB(std::sin(B)), cosB(std::cos(B));

                // p1
                const float_type I(1), S(a);
                const float_type dy((inv ? -1 : 1) * S / (I * a_0) * dB);
                const float_type dx(dy * (p1.y - l * sinB) / (l * cosB - r + F));
                const float_vec p1_new(p1.x + dx, p1.y + dy);

                // p2
                const float_type lump((-R * K_out.y + l * sinB - p1.y) / (R * K_out.x - l * cosB + p1.x));
                const float_type dl(l * (lump * cosB + sinB) / (cosB - lump * sinB) * dB);
                const float_type l_new(l + dl);
                const float_vec p2_new(l_new * float_vec(std::cos(B_new), std::sin(B_new)));

                // directions
                float_vec K_int(p2_new - p1_new);
                const float_type R_new(length(K_int));
                K_int /= R_new;
                auto cone(calcCone(p2_new, tpts));
                a = cone.first;
                K_out = cone.second;
                if (HitInfo hh; traceRay(Ray(p2, K_out), hh, t)) r = hh.l;
                else std::cout << "error: no hit" << std::endl;

                // normals
                const float_vec n1_new(normalize(K_int - K_in));
                const float_vec n2_new(normalize(K_out - K_int));

                // store coordinates and normals
                m1a.addSegment(p1, p1_new, n1, n1_new);
                m1b.addSegment(float_vec(p1.x, -p1.y), float_vec(p1_new.x, -p1_new.y), float_vec(n1.x, -n1.y), float_vec(n1_new.x, -n1_new.y));
                m2a.addSegment(p2, p2_new, n2, n2_new);
                m2b.addSegment(float_vec(p2.x, -p2.y), float_vec(p2_new.x, -p2_new.y), float_vec(n2.x, -n2.y), float_vec(n2_new.x, -n2_new.y));

                // update
                B = B_new;
                p1 = p1_new;
                p2 = p2_new;
                n1 = n1_new;
                n2 = n2_new;
                R = R_new;
                l = l_new;
            }
        }

        // build barrier
        buildBarrier();
    }



    void buildBarrier() {
        const float_type x_max = 1.25 * std::max(std::abs(m1a.segments.front().p1.x), std::abs(m2a.segments.front().p1.x));
        const float_type y_max = 1.25 * std::max(std::abs(m1a.segments.back().p2.y), std::abs(m2a.segments.back().p2.y));
        const float_vec bottomRight(x_max, -y_max), topRight(x_max, y_max), topLeft(-x_max, y_max), bottomLeft(-x_max, -y_max);
        barrier.addSegment(bottomRight, topRight, float_vec(-1, 0), float_vec(-1, 0));
        barrier.addSegment(topRight, topLeft, float_vec(0, -1), float_vec(0, -1));
        barrier.addSegment(topLeft, bottomLeft, float_vec(1, 0), float_vec(1, 0));
        barrier.addSegment(bottomLeft, bottomRight, float_vec(0, 1), float_vec(0, 1));
    }

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override {
        if (m1a.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_1a;
            return true;
        }
        if (m1b.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_1b;
            return true;
        }
        if (m2a.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_2a;
            return true;
        }
        if (m2b.intersect(ray, hitInfo)) {
            hitInfo.t = Type::MIRROR_2b;
            return true;
        }
        if (barrier.intersect(ray, hitInfo)) {
            hitInfo.t = Type::BARRIER;
            return true;
        }
        return false;
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

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
    }

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override {
        return std::pair(std::vector<Ray>(), std::vector<Ray>());
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



// system
struct Design {

    Geometry* source;
    Geometry* target;
    std::vector<Geometry*> geometries;

    void addGeometry(Geometry* const g) {
        geometries.push_back(g);
        if (g->type == Type::SOURCE) source = g;
        else if (g->type == Type::TARGET) target = g;
    }

    bool intersect(const Ray& ray, HitInfo& minHitInfo) const {
        bool hit = false;
        HitInfo tempMinHitInfo;
        minHitInfo.l = std::numeric_limits<float_type>::max();
        for (const auto& geometry : geometries) {
            if ((*geometry).intersect(ray, tempMinHitInfo)) {
                if (tempMinHitInfo.l < minHitInfo.l) {
                    hit = true;
                    minHitInfo = tempMinHitInfo;
                }
            }
        }
        return hit;
    }

    Path traceRay(const Ray& ray) const {
        Path path;
        path.addVertex(ray.o);
        Ray r = ray;
        for (int i = 0; i < 3; ++i) {
            if (HitInfo h; intersect(r, h)) {
                path.addVertex(h.p);
                if (h.t == Type::TARGET) break;
                float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
            }
            else break;
        }
        return path;
    }

    std::vector<Path> rayTrace(const std::vector<Ray>& rays) const {
        std::vector<std::vector<Path>> threadVectors(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (const auto& ray : rays) {
            int threadID = omp_get_thread_num();
            threadVectors[threadID].push_back(traceRay(ray));
        }
        std::vector<Path> paths;
        for (const auto& v : threadVectors) paths.insert(paths.end(), v.begin(), v.end());
        return paths;
    }

    void traceMeanRays(const std::string& filePath, const int& numRays) {
        std::vector<Ray> rays;
        for (int i = 0; i < numRays; ++i) rays.push_back(source->sampleMeanRay());
        std::vector<Path> paths(rayTrace(rays));
        std::ofstream file;
        file.open(filePath + "mean.csv");
        for (auto path : paths) path.writePath(file);
        file.close();
    }

    void traceDiffuseRays(const std::string& filePath, const int& numRays) const {
        std::vector<Ray> rays;
        for (int i = 0; i < numRays; ++i) rays.push_back(source->sampleDiffuseRay().first);
        std::vector<Path> paths(rayTrace(rays));
        std::ofstream file;
        file.open(filePath + "diffuse.csv");
        for (auto path : paths) path.writePath(file);
        file.close();
    }

    void traceExtremeDiffuseRays(const std::string& filePath) const {
        auto extremeRays = source->generateExtremeDiffuseRays();
        auto paths = std::pair(rayTrace(extremeRays.first), rayTrace(extremeRays.second));
        std::ofstream file;
        file.open(filePath + "extrdiff1.csv");
        for (auto path : paths.first) path.writeFinalSegment(file);
        file.close();
        file.open(filePath + "extrdiff2.csv");
        for (auto path : paths.second) path.writeFinalSegment(file);
        file.close();
    }

    void traceExtremeInfiniteRays(const std::string& filePath, const int& numRays) {
        auto extremeRays = source->generateExtremeInfiniteRays(numRays);
        auto paths = std::pair(rayTrace(extremeRays.first), rayTrace(extremeRays.second));
        std::ofstream file;
        file.open(filePath + "extrinf1.csv");
        for (auto path : paths.first) path.writeFinalSegment(file);
        file.close();
        file.open(filePath + "extrinf2.csv");
        for (auto path : paths.second) path.writeFinalSegment(file);
        file.close();
    }



    void tracePhaseSpace(const std::string& filePath, const int& numRays) const {

        std::vector<std::vector<std::vector<float_type>>> threadVectorsTarget(omp_get_max_threads());
        std::vector<std::vector<std::vector<float_type>>> threadVectorsMirror1(omp_get_max_threads());
        std::vector<std::vector<std::vector<float_type>>> threadVectorsMirror2(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            auto sample = source->sampleDiffuseRay();
            Ray r = sample.first;
            const float_type source_angle(sample.second);
            const float_vec source_pos(r.o - source->getCentre());

            // for 3 path segments
            for (int j = 0; j < 3; ++j) {

                // check intersection
                if (HitInfo h; intersect(r, h)) {

                    // mirror 1
                    if (j == 0 && h.t != Type::MIRROR_1a && h.t != Type::MIRROR_1b) break;
                    if (j == 0 && (h.t == Type::MIRROR_1a || h.t == Type::MIRROR_1b)) {
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, std::abs(source_pos.y)};
                        threadVectorsMirror1[threadID].emplace_back(std::move(temp));
                    }

                    // mirror 2
                    if (j == 1 && h.t != Type::MIRROR_2a && h.t != Type::MIRROR_2b) break;
                    if (j == 1 && (h.t == Type::MIRROR_2a || h.t == Type::MIRROR_2b)) {
                        std::vector<float_type> temp = {cross(-r.d, h.n), h.ml, source_angle, std::abs(source_pos.y)};
                        threadVectorsMirror2[threadID].emplace_back(std::move(temp));
                    }

                    // target
                    if (j == 2 && h.t != Type::TARGET) break;
                    if (j == 2 && h.t == Type::TARGET) {
                        if (target->shape == Shape::FLAT) {
                            std::vector<float_type> temp = {cross(-r.d, h.n), h.p.y, source_angle, std::abs(source_pos.y)};
                            threadVectorsTarget[threadID].emplace_back(std::move(temp));
                        }
                        else if (target->shape == Shape::CYLINDRICAL) {
                            float_type source_theta = std::atan2(source_pos.y, -source_pos.x);
                            if (source_theta < 0) source_theta += 2 * PI;
                            float_type target_theta = std::atan2(h.p.y, h.p.x);
                            if (target_theta < 0) target_theta += 2 * PI;
                            std::vector<float_type> temp = {cross(-r.d, h.n), target_theta * RAD_TO_DEG, source_angle, source_theta * RAD_TO_DEG};
                            threadVectorsTarget[threadID].emplace_back(std::move(temp));
                        }
                        break;
                    }

                    // reflect
                    float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                    r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
                }
                else break;
            }
        }

        std::ofstream file;

        // target
        std::vector<std::vector<float_type>> dataTarget;
        for (const auto& v : threadVectorsTarget) dataTarget.insert(dataTarget.end(), v.begin(), v.end());
        file.open(filePath + "phasetarget.csv");
        for (const auto& p : dataTarget) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
        file.close();

        // mirror 1
        std::vector<std::vector<float_type>> dataMirror1;
        for (const auto& v : threadVectorsMirror1) dataMirror1.insert(dataMirror1.end(), v.begin(), v.end());
        file.open(filePath + "phasemirror1.csv");
        for (const auto& p : dataMirror1) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
        file.close();

        // mirror 2
        std::vector<std::vector<float_type>> dataMirror2;
        for (const auto& v : threadVectorsMirror2) dataMirror2.insert(dataMirror2.end(), v.begin(), v.end());
        file.open(filePath + "phasemirror2.csv");
        for (const auto& p : dataMirror2) file << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "\n";
        file.close();
    }



    float_type traceHitData(const int& numRays) const {
        std::vector<int> threadCounts(omp_get_max_threads());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            Ray r = source->sampleDiffuseRay().first;
            for (int j = 0; j < 3; ++j) {
                if (HitInfo h; intersect(r, h)) {
                    if (h.t == Type::TARGET) {
                        threadCounts[threadID] += 1;
                        break;
                    }
                    if (h.t == Type::BARRIER) break;
                    // if (j == 0 && h.t != Type::MIRROR_1a && h.t != Type::MIRROR_1b) break;
                    // if (j == 1 && h.t != Type::MIRROR_2a && h.t != Type::MIRROR_2b) break;
                    // if (j == 2 && h.t != Type::TARGET) break;
                    // if (j == 2 && h.t == Type::TARGET) {
                    //     threadCounts[threadID] += 1;
                    //     break;
                    // }
                    float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                    r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
                }
                else break;
            }
        }
        int finalCount = 0;
        for (const auto& count : threadCounts) finalCount += count;
        return static_cast<float_type>(finalCount) / static_cast<float_type>(numRays);
    }

    std::vector<float_type> traceDetailedHitData(const int& numRays) const {
        std::vector<std::vector<int>> threadCounts(omp_get_max_threads(), std::vector<int>(4, 0));
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < numRays; ++i) {
            int threadID = omp_get_thread_num();
            Ray r = source->sampleDiffuseRay().first;
            for (int j = 0; j < 3; ++j) {
                if (HitInfo h; intersect(r, h)) {
                    if (h.t == Type::TARGET) {
                        if (j == 0) threadCounts[threadID][0] += 1;
                        else if (j == 1) threadCounts[threadID][1] += 1;
                        else if (j == 2) threadCounts[threadID][2] += 1;
                        break;
                    }
                    else if (h.t == Type::BARRIER) {
                        threadCounts[threadID][3] += 1;
                        break;
                    }
                    float_vec refl = normalize(r.d - 2 * dot(r.d, h.n) * h.n);
                    r = Ray(h.p + (DOINKING ? DOINK * refl : float_vec(0, 0)), refl);
                }
                else break;
            }
        }
        std::vector<int> finalCounts(4, 0);
        for (const auto& counts : threadCounts) {
            finalCounts[0] += counts[0];
            finalCounts[1] += counts[1];
            finalCounts[2] += counts[2];
            finalCounts[3] += counts[3];
        }
        std::vector<float_type> finalMetrics(4);
        finalMetrics[0] = static_cast<float_type>(finalCounts[0]) / static_cast<float_type>(numRays);
        finalMetrics[1] = static_cast<float_type>(finalCounts[1]) / static_cast<float_type>(numRays);
        finalMetrics[2] = static_cast<float_type>(finalCounts[2]) / static_cast<float_type>(numRays);
        finalMetrics[3] = static_cast<float_type>(finalCounts[3]) / static_cast<float_type>(numRays);
        return finalMetrics;
    }

};



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



// main
int main(int argc, char* argv[]) {



    // finite system
    if (BUILD_FINITE_SYSTEM) {

        // output
        std::string outputDataPath = "datafin/";
        if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
        std::ofstream file;

        // input
        const bool inv(false);
        const float_type f1(1), L(4), f2(1), da(0.0001), a_max(90 * DEG_TO_RAD), radius(f1 / 10);
        file.open(outputDataPath + "input.csv");
        file << FLAT_SOURCE << "," << FLAT_TARGET << "," << inv << "," << radius << "," << f1 << "," << L << "," << f2 << "," << da << "," << a_max * RAD_TO_DEG << "\n";
        file.close();

        // source
        LineSegment flatSource(Type::SOURCE, float_vec(-L, -radius), float_vec(-L, radius), float_vec(-1, 0), float_vec(-1, 0));
        Circle cylindricalSource(Type::SOURCE, float_vec(-L, 0), radius);

        // target
        LineSegment flatTarget(Type::TARGET, float_vec(0, -radius), float_vec(0, radius), float_vec(1, 0), float_vec(1, 0));
        Circle cylindricalTarget(Type::TARGET, float_vec(0, 0), radius);

        // two-mirror concentrator
        TwoMirrorConcentrator tmc;
        tmc.buildFin(inv, f1, L, f2, da, a_max);
        tmc.writeTwoMirrorConcentrator(outputDataPath);

        // design
        Design design;
        if (FLAT_SOURCE) design.addGeometry(&flatSource);
        else design.addGeometry(&cylindricalSource);
        design.addGeometry(&tmc);

        // extreme rays
        design.traceExtremeDiffuseRays(outputDataPath);

        // add target
        if (FLAT_TARGET) design.addGeometry(&flatTarget);
        else design.addGeometry(&cylindricalTarget);

        // phase space
        design.tracePhaseSpace(outputDataPath, 10000);

        /*

        // hit report
        const int numberOfDesigns = 20;
        std::vector<Design> designs(numberOfDesigns);
        std::vector<LineSegment> flatTargets;
        std::vector<Circle> cylindricalTargets;
        const float_type increment = 2 * radius / numberOfDesigns;
        for (int i = 0; i < numberOfDesigns; ++i) {
            if (FLAT_TARGET) flatTargets.emplace_back(Type::TARGET, float_vec(0, -(i + 1) * increment), float_vec(0, (i + 1) * increment), float_vec(1, 0), float_vec(1, 0));
            else cylindricalTargets.emplace_back(Type::TARGET, float_vec(0, 0), (i + 1) * increment);
        }
        for (int i = 0; i < numberOfDesigns; ++i) {
            if (FLAT_SOURCE) designs[i].addGeometry(&flatSource);
            else designs[i].addGeometry(&cylindricalSource);
            if (FLAT_TARGET) designs[i].addGeometry(&flatTargets[i]);
            else designs[i].addGeometry(&cylindricalTargets[i]);
            designs[i].addGeometry(&tmc);
        }
        std::vector<float_vec> hitData;
        const int numberOfTrials = 3;
        const int numberOfRays = 10000;
        for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitData(numberOfRays));
        file.open(outputDataPath + "hitdata.csv");
        for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
        file.close();

        // std::vector<std::pair<float_type, std::vector<float_type>>> hitData;
        // const int numberOfTrials = 3;
        // const int numberOfRays = 10000;
        // for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceDetailedHitData(numberOfRays));
        // file.open(outputDataPath + "hitdata.csv");
        // for (const auto& p : hitData) file << p.first << "," << p.second[0] << "," << p.second[1] << "," << p.second[2] << "," << p.second[3] << "\n";
        // file.close();

        */

    }



    // infinite system
    if (BUILD_INFINITE_SYSTEM) {

        // output
        std::string outputDataPath = "datainf/";
        if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
        std::ofstream file;

        // input
        const bool inv(true);
        const float_type L(3), f(0.5), dB(0.000001), B_max(85 * DEG_TO_RAD);
        const float_vec K_in(-1, 0);
        file.open(outputDataPath + "input.csv");
        file << FLAT_TARGET << "," << inv << "," << L << "," << f << "," << dB << "," << B_max * RAD_TO_DEG << "\n";
        file.close();

        // two-mirror concentrator
        TwoMirrorConcentrator tmc;
        tmc.buildInf(inv, L, f, K_in, dB, B_max);
        tmc.writeTwoMirrorConcentrator(outputDataPath);

        // design
        Design d;
        d.addGeometry(&tmc);

        // infinite source
        LineSegment source(Type::SOURCE, float_vec(0, -2), float_vec(0, 2), K_in, K_in);
        // LineSegment source(Type::SOURCE, float_vec(-1, -2), float_vec(-1, 2), K_in, K_in);
        d.addGeometry(&source);

        // trace extreme rays
        d.traceExtremeInfiniteRays(outputDataPath, 50);

        // flat target
        LineSegment target(Type::TARGET, float_vec(0, -EPSILON_OVER_TWO), float_vec(0, EPSILON_OVER_TWO), float_vec(1, 0), float_vec(1, 0));
        d.addGeometry(&target);
    }



    // infinite system w/ arbitrary target
    if (BUILD_INFINITE_ARBITRARY_SYSTEM) {

        // output
        std::string outputDataPath = "datainfarb/";
        if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
        std::ofstream file;

        // input
        const bool inv(true);
        const float_type L(3), f(0.5), dB(0.000001), B_max(85 * DEG_TO_RAD);
        const float_vec K_in(-1, 0);
        file.open(outputDataPath + "input.csv");
        file << inv << "," << L << "," << f << "," << dB << "," << B_max * RAD_TO_DEG << "\n";
        file.close();

        // two-mirror concentrator
        TwoMirrorConcentrator tmc;
        tmc.buildInfArb(inv, L, f, K_in, dB, B_max);
        tmc.writeTwoMirrorConcentrator(outputDataPath);

        // design
        Design d;
        d.addGeometry(&tmc);

        // infinite source
        LineSegment source(Type::SOURCE, float_vec(0, -2), float_vec(0, 2), K_in, K_in);
        // LineSegment source(Type::SOURCE, float_vec(-1, -2), float_vec(-1, 2), K_in, K_in);
        d.addGeometry(&source);

        // trace extreme rays
        d.traceExtremeInfiniteRays(outputDataPath, 50);

        // flat target
        LineSegment target(Type::TARGET, float_vec(0, -EPSILON_OVER_TWO), float_vec(0, EPSILON_OVER_TWO), float_vec(1, 0), float_vec(1, 0));
        d.addGeometry(&target);

        /*

        // triangular target
        // const float_vec tp1(1, 1), tp2(-2, 0), tp3(0, -1);
        // const float_vec tt1(tp1 - tp2), tt2(tp2 - tp3), tt3(tp3 - tp1);
        // const float_vec tn1(-tt1.y, tt1.x), tn2(-tt2.y, tt2.x), tn3(-tt3.y, tt3.x);
        // LineSegment t1(Type::TARGET, tp1, tp2, tn1, tn1);
        // LineSegment t2(Type::TARGET, tp2, tp3, tn2, tn2);
        // LineSegment t3(Type::TARGET, tp3, tp1, tn3, tn3);
        // d.addGeometry(&t1);
        // d.addGeometry(&t2);
        // d.addGeometry(&t3);

        */
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

/*

// ray trace
std::vector<Ray> rays;
if (FLAT_SOURCE) {
    for (int i = 0; i < 100; ++i) {
        // rays.push_back(flatSource.sampleMeanRay());
        rays.push_back(flatSource.sampleDiffuseRay());
    }
}
if (CYLINDRICAL_SOURCE) {
    for (int i = 0; i < 100; ++i) {
        // rays.push_back(cylindricalSource.sampleMeanRay());
        rays.push_back(cylindricalSource.sampleDiffuseRay());
    }
}
std::vector<Path> paths = design.rayTrace(rays);

// sampled rays
file.open("data_2mc_fin_paths.csv");
for (auto path : paths) {
    path.writePath(file);
}
file.close();

*/

/*

// extreme rays
tmc.buildBarrier();
if (FLAT_SOURCE) {
    auto extrRays = flatSource.generateExtremeRays();
    auto p1Paths = design.rayTrace(extrRays.first);
    auto p2Paths = design.rayTrace(extrRays.second);
    file.open(outputDataPath + "extreme1.csv");
    for (auto path : p1Paths) path.writeFinalSegment(file);
    file.close();
    file.open(outputDataPath + "extreme2.csv");
    for (auto path : p2Paths) path.writeFinalSegment(file);
    file.close();
}
if (CYLINDRICAL_SOURCE) {
    auto extrRays = cylindricalSource.generateExtremeRays();
    auto p1Paths = design.rayTrace(extrRays.first);
    auto p2Paths = design.rayTrace(extrRays.second);
    file.open(outputDataPath + "extreme1.csv");
    for (auto path : p1Paths) path.writeFinalSegment(file);
    file.close();
    file.open(outputDataPath + "extreme2.csv");
    for (auto path : p2Paths) path.writeFinalSegment(file);
    file.close();
}

// phase space
if (FLAT_SOURCE) design.addGeometry(&flatSource);
else if (CYLINDRICAL_SOURCE) design.addGeometry(&cylindricalSource);
if (FLAT_TARGET) design.addGeometry(&flatTarget);
else if (CYLINDRICAL_TARGET) design.addGeometry(&cylindricalTarget);
auto phaseSpaceData = design.tracePhaseSpace(10000);
file.open(outputDataPath + "phase.csv");
for (auto data : phaseSpaceData) file << data[0] << "," << data[1] << "," << data[2] << "\n";
file.close();

*/
