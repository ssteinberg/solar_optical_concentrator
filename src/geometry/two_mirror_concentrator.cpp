#include "geometry/two_mirror_concentrator.h"
#include "pcg32.h"

#include <fstream>
#include <iostream>

void TwoMirrorConcentrator::reset() {
    m1a.reset();
    m1b.reset();
    m2a.reset();
    m2b.reset();
    barrier.reset();
}

float_vec TwoMirrorConcentrator::getCentre() const {
    return float_vec(0, 0);
}

float_type TwoMirrorConcentrator::getLength() const {
    return m1a.getLength() + m1b.getLength() + m2a.getLength() + m2b.getLength();
}

void TwoMirrorConcentrator::buildFin(const bool &inv, const float_type &f1, const float_type &L, const float_type &f2,
    const float_type &da, const float_type &a_max, const float_type &w) {

    // reset
    reset();

    // initial conditions
    const auto& Sa = [=](const float_type& alpha) {
        if (FLAT_SOURCE) return std::cos(alpha);
        else return float_type(1);
    };
    const auto& SB = [=](const float_type& beta) {
        if constexpr (FLAT_SOURCE && TARGET_SHAPE == Shape::FLAT) return std::cos(beta) * w;
        if constexpr (FLAT_SOURCE && TARGET_SHAPE == Shape::CYLINDRICAL) return ONE_OVER_PI;
        if constexpr (!FLAT_SOURCE && TARGET_SHAPE == Shape::FLAT) return PI * std::cos(beta);
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

void TwoMirrorConcentrator::buildInf(const bool &inv, const float_type &L, const float_type &f, const float_vec &K_in,
    const float_type &dB, const float_type &B_max) {

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
        const float_type I(1), S(getAngularIntensityDistribution(B));
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

bool TwoMirrorConcentrator::calcCone(const Polygon &p, const float_vec &apex, HitInfo &h) {

    // find u
    float_type sina_max(0);
    float_vec u(0, 0);
    for (const auto& v1 : p.vertices) {
        for (const auto& v2 : p.vertices) {
            const float_vec l1(normalize(v1 - apex)), l2(normalize(v2 - apex));
            const float_type sina(std::abs(cross(l1, l2)));
            if (sina > sina_max) {
                sina_max = sina;
                u = normalize(l1 + l2);
            }
        }
    }

    // find rtraced
    bool hit = false;
    HitInfo temp;
    h.l = std::numeric_limits<float_type>::max();
    if(p.intersect(Ray(apex, u), temp)) {
        if (temp.l < h.l) {
            hit = true;
            h = temp;
        }
    }
    if (hit) h.rtraced = h.l;
    else return false;

    // find rmin, rmax, rmean
    float_type rmin(std::numeric_limits<float_type>::max()), rmax(std::numeric_limits<float_type>::min());
    float_type rmin2(std::numeric_limits<float_type>::max()), rmax2(std::numeric_limits<float_type>::min());
    for (const auto& v : p.vertices) {
        const float_vec seg(v - apex);
        const float_type proj(dot(seg, u));
        if (proj > rmax) rmax = proj;
        if (proj < rmin) rmin = proj;
        const float_type len(length(seg));
        if (len > rmax2) rmax2 = len;
        if (len < rmin2) rmin2 = len;
    }
    h.rmin = rmin;
    h.rmax = rmax;
    h.rmean = 0.5 * (rmin + rmax);
    h.rmin2 = rmin2;
    h.rmax2 = rmax2;
    h.rmean2 = 0.5 * (rmin2 + rmax2);

    // return
    h.u = u;
    h.a = std::asin(sina_max);
    h.sina = sina_max;
    h.tana = 2 * std::tan(h.a / 2);
    return true;
}

void TwoMirrorConcentrator::buildInfArb(const Polygon &p, const bool &inv, const float_type &L, const float_type &f,
    const float_vec &K_in, const float_type &dB, const float_type &B_max, const int &i, const int &j) {

    // reset
    reset();

    // initial conditions
    float_vec p1(-L, 0), p2(f, 0), n1(1, 0);
    float_type R(L + f), l(f);

    // calculate cone to initialize
    if (HitInfo h; calcCone(p, p2, h)) {
        float_vec K_out(h.u), n2(K_out - normalize(p2 - p1));
        float_type r_0, r, a_0, a;
        if (i == 0) r_0 = r = h.rtraced;
        else if (i == 1) r_0 = r = h.rmin;
        else if (i == 2) r_0 = r = h.rmax;
        else if (i == 3) r_0 = r = h.rmean;
        else if (i == 4) r_0 = r = h.rmin2;
        else if (i == 5) r_0 = r = h.rmax2;
        else r_0 = r = h.rmean2;
        if (j == 0) a_0 = a = h.a;
        else if (j == 1) a_0 = a = h.sina;
        else a_0 = a = r * h.tana;
        const float_type F(2 * L + f + r_0);

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
            if (HitInfo hh; calcCone(p, p2, hh)) {
                K_out = hh.u;
                if (i == 0) r = hh.rtraced;
                else if (i == 1) r = hh.rmin;
                else if (i == 2) r = hh.rmax;
                else if (i == 3) r = hh.rmean;
                else if (i == 4) r = hh.rmin2;
                else if (i == 5) r = hh.rmax2;
                else r = hh.rmean2;
                if (j == 0) a = hh.a;
                else if (j == 1) a = hh.sina;
                else a = r * hh.tana;
            }
            else {
                std::cout << "ERROR: No hit, breaking loop." << std::endl;
                break;
            }

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

void TwoMirrorConcentrator::buildBarrier() {
    const float_type x_max = 1.25 * std::max(std::abs(m1a.segments.front().p1.x), std::abs(m2a.segments.front().p1.x));
    const float_type y_max = 1.25 * std::max(std::abs(m1a.segments.back().p2.y), std::abs(m2a.segments.back().p2.y));
    const float_vec bottomRight(x_max, -y_max), topRight(x_max, y_max), topLeft(-x_max, y_max), bottomLeft(-x_max, -y_max);
    barrier.addSegment(bottomRight, topRight, float_vec(-1, 0), float_vec(-1, 0));
    barrier.addSegment(topRight, topLeft, float_vec(0, -1), float_vec(0, -1));
    barrier.addSegment(topLeft, bottomLeft, float_vec(1, 0), float_vec(1, 0));
    barrier.addSegment(bottomLeft, bottomRight, float_vec(0, 1), float_vec(0, 1));
}

bool TwoMirrorConcentrator::intersect(const Ray &ray, HitInfo &hitInfo) const {
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

Ray TwoMirrorConcentrator::sampleMeanRay() const {
    const float_type Randy = PCG32::rand();
    if (Randy < 0.25) return m1a.sampleMeanRay();
    if (0.25 <= Randy && Randy < 0.5) return m1b.sampleMeanRay();
    if (0.5 <= Randy && Randy < 0.75) return m2a.sampleMeanRay();
    return m2b.sampleMeanRay();
}

std::pair<Ray, float_type> TwoMirrorConcentrator::sampleDiffuseRay() const {
    const float_type Randolf = PCG32::rand();
    if (Randolf < 0.25) return m1a.sampleDiffuseRay();
    if (0.25 <= Randolf && Randolf < 0.5) return m1b.sampleDiffuseRay();
    if (0.5 <= Randolf && Randolf < 0.75) return m2a.sampleDiffuseRay();
    return m2b.sampleDiffuseRay();
}

std::pair<std::vector<Ray>, std::vector<Ray>> TwoMirrorConcentrator::generateExtremeDiffuseRays() const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::pair<std::vector<Ray>, std::vector<Ray>> TwoMirrorConcentrator::generateExtremeInfiniteRays(
    const int &numRays) const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::vector<Ray> TwoMirrorConcentrator::generateFinalPlotRays() const {
    return std::vector<Ray>();
}

void TwoMirrorConcentrator::writeTwoMirrorConcentrator(const std::string &filePath) const {
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

float_type TwoMirrorConcentrator::getAngularIntensityDistribution(const float_type beta) {
    float_type angularIntensityDistribution { 0 };
    switch (TARGET_SHAPE) {
        case Shape::FLAT:
            angularIntensityDistribution = std::cos(beta);
            break;
        case Shape::CYLINDRICAL:
            angularIntensityDistribution = 1;
            break;
        case Shape::ELLIPTICAL:
            angularIntensityDistribution = std::sqrt(std::pow(ELLIPTICAL_TARGET_X_RADIUS, 2) * std::pow(std::sin(beta), 2) + std::pow(ELLIPTICAL_TARGET_Y_RADIUS, 2) * std::pow(std::cos(beta), 2)) / std::max(ELLIPTICAL_TARGET_X_RADIUS, ELLIPTICAL_TARGET_Y_RADIUS);
            break;
        default:
            throw std::invalid_argument("Target type is not supported.");
    }

    if constexpr (SHOULD_WEIGHT_ANGULAR_INTENSITY_DISTRIBUTION) {
        angularIntensityDistribution *= getOmega(beta);
    }

    return angularIntensityDistribution;
}

// Weighting for the desired angular intensity distribution
float_type TwoMirrorConcentrator::getOmega(const float_type beta) {
    return std::exp(-(beta * beta)/(2 * ANGULAR_INTENSITY_DISTRIBUTION_WEIGHTING_C * ANGULAR_INTENSITY_DISTRIBUTION_WEIGHTING_C));
}
