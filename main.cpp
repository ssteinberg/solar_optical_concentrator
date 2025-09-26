////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//  compile: "compiler -std=c++17 main.cpp", compiler = "g++-15", "g++", "clang++"
//      run: "./a.out"
//     plot: run results.ipynb

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// standard library
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <cmath>

// linear algebra
#include "linalg.h"
using namespace linalg::aliases;

// floating-point precision
using float_type = double;
constexpr int NUMBER_OF_DIMENSIONS = 2;
typedef linalg::vec<float_type, NUMBER_OF_DIMENSIONS> float_vec;

// global constants
constexpr float_type PI = 3.14159265358979;
constexpr float_type PI_OVER_TWO = PI / 2;
constexpr float_type DEG_TO_RAD = PI / 180;
constexpr float_type RAD_TO_DEG = 180 / PI;
constexpr float_type EPSILON = 0.01;
constexpr float_type EPSILON_OVER_TWO = EPSILON / 2;

// choose what to build
constexpr bool INFINITE_SYSTEM = false;
constexpr bool FINITE_SYSTEM = true;
constexpr bool FLAT_TARGET = false;
constexpr bool CYLINDRICAL_TARGET = false;
constexpr bool PARABOLIC_MIRROR = false;
constexpr bool SCALE_INVARIANT = false;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// ray
struct Ray {
    float_vec o, d;
    Ray(const float_vec& o, const float_vec& d) : o(o), d(d) {}
};

// hit info
struct HitInfo {
    float_vec P, N;
};

// target shape
enum struct TargetShape {
    FLAT,
    CYLINDRICAL
};

// mirror
struct Mirror {

    // coordinates and normals of the line segments of the mirror
    std::vector<std::pair<float_vec, float_vec>> segments;
    std::vector<std::pair<float_vec, float_vec>> normals;

    // add coordinates and normals of a line segment to the mirror
    void addSegment(const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2) {
        segments.emplace_back(p1, p2);
        normals.emplace_back(n1, n2);
    }

    // check for ray-mirror intersection
    bool intersect(HitInfo& hitInfo, const Ray& ray) {

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
            if (-0.00001 <= t && t <= 1.00001) {
                hitInfo.P = segments[i].first + t * segmentDirection;
                hitInfo.N = normalize(t * normals[i].second + (1 - t) * normals[i].first);
                return true;
            }
        }

        // no intersection
        return false;
    }

};

// two-mirror concentrator
struct TwoMirrorConcentrator {

    // fields
    Mirror M1, M2;

    // build infinite system
    void buildInfinite(const TargetShape& ts, const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max) {

        // incoming intensity
        const auto& I = [=](const float_type alpha) {
            if (ts == TargetShape::FLAT) return std::cos(alpha);
            else if (ts == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

        // outgoing intensity
        const auto& S = [=](const float_type beta) {
            if (ts == TargetShape::FLAT) return std::cos(beta);
            else if (ts == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

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
            const float_type I = 1; // depends on I = const.
            const float_type dy = (inv ? -1 : 1) * S(B) / I * dB;
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
            M1.addSegment(p, p_new, n, n_new);
            M2.addSegment(pp, pp_new, np, np_new);

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

    // build finite system
    void buildFinite(const TargetShape& ts, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max) {

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
            const float_type dB = da;
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
            M1.addSegment(p1, p1_new, n1, n1_new);
            M2.addSegment(p2, p2_new, n2, n2_new);

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

    // trace the path of a ray in the system
    void traceRay(std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>>& rayPaths, const Ray& ray1) {

        // check intersection with mirror 1
        if (HitInfo hit1; M1.intersect(hit1, ray1)) {

            // check intersection with mirror 2
            const Ray ray2(hit1.P, normalize(ray1.d - 2 * dot(hit1.N, ray1.d) * hit1.N));
            if (HitInfo hit2; M2.intersect(hit2, ray2)) {

                // compute point through target and store path of ray
                const float_type len = 6;
                rayPaths.emplace_back(ray1.o, hit1.P, hit2.P, hit2.P + len * normalize(ray2.d - 2 * dot(hit2.N, ray2.d) * hit2.N));
            }
        }
    }

    // trace multiple rays across y direction
    void traceRays(std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>>& rayPaths, const float_vec& rayDirection,
        const float_type& h_min, const float_type& h_inc, const float_type& h_max, const std::string& fileName) {

        // trace
        float_type h = h_min;
        while (h < h_max) {
            const Ray ray(float_vec(0, h), rayDirection);
            traceRay(rayPaths, ray);
            h += h_inc;
        }

        // write
        std::ofstream file;
        file.open(fileName);
        for (auto rayPath : rayPaths) {
            file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
            file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
            file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
            file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
        }
        file.close();
    }

};

// parabolic mirror
struct ParabolicMirror {

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
        if (HitInfo hit; M.intersect(hit, ray)) {
            rayPaths.emplace_back(ray.o, hit.P, hit.P + normalize(ray.d - 2 * dot(hit.N, ray.d) * hit.N));
        }
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// build infinite two-mirror concentrator
TwoMirrorConcentrator buildInf2MC(const TargetShape& ts, const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max,
    const std::string& fileName) {

    // build
    TwoMirrorConcentrator inf2MC;
    inf2MC.buildInfinite(ts, inv, L, f, K_in, dB, B_max);

    // write
    std::ofstream file;
    file.open(fileName);
    for (int i = 0, n = static_cast<int>(inf2MC.M1.segments.size()); i < n; ++i) {
        file << inf2MC.M1.segments[i].first.x << "," << inf2MC.M1.segments[i].first.y << ",";
        file << inf2MC.M2.segments[i].first.x << "," << inf2MC.M2.segments[i].first.y << "\n";
    }
    file << inf2MC.M1.segments.back().second.x << "," << inf2MC.M1.segments.back().second.y << ",";
    file << inf2MC.M2.segments.back().second.x << "," << inf2MC.M2.segments.back().second.y << "\n";
    file.close();

    // return
    return inf2MC;
}

// build finite two-mirror concentrator
TwoMirrorConcentrator buildFin2MC(const TargetShape& ts, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max,
    const std::string& fileName) {
    
    // build
    TwoMirrorConcentrator fin2MC;
    fin2MC.buildFinite(ts, f1, L, f2, da, a_max);

    // write
    std::ofstream file;
    file.open(fileName);
    for (int i = 0, n = static_cast<int>(fin2MC.M1.segments.size()); i < n; ++i) {
        file << fin2MC.M1.segments[i].first.x << "," << fin2MC.M1.segments[i].first.y << ",";
        file << fin2MC.M2.segments[i].first.x << "," << fin2MC.M2.segments[i].first.y << "\n";
    }
    file << fin2MC.M1.segments.back().second.x << "," << fin2MC.M1.segments.back().second.y << ",";
    file << fin2MC.M2.segments.back().second.x << "," << fin2MC.M2.segments.back().second.y << "\n";
    file.close();

    // return
    return fin2MC;
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
    for (auto segment : PM.M.segments) file << segment.first.x << "," << segment.first.y << "\n";
    file.close();

    // return
    return PM;
}



// main
int main(int argc, char* argv[]) {

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
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> meanPaths;
            infFlat2MC.traceRays(meanPaths, K_in, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_mean.csv");
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> extrPosPaths;
            float_vec extrPosDir(-std::cos(EPSILON_OVER_TWO), std::sin(EPSILON_OVER_TWO));
            infFlat2MC.traceRays(extrPosPaths, extrPosDir, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_extr_pos.csv");
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> extrNegPaths;
            float_vec extrNegDir(-std::cos(EPSILON_OVER_TWO), -std::sin(EPSILON_OVER_TWO));
            infFlat2MC.traceRays(extrNegPaths, extrNegDir, -infFlat2MC.M1.segments.back().second.y, h_inc, infFlat2MC.M1.segments.back().second.y, "data_2mc_flat_extr_neg.csv");
        }
        if (CYLINDRICAL_TARGET) {
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> meanPaths;
            infCyl2MC.traceRays(meanPaths, K_in, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_mean.csv");
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> extrPosPaths;
            float_vec extrPosDir(-std::cos(EPSILON_OVER_TWO), std::sin(EPSILON_OVER_TWO));
            infCyl2MC.traceRays(extrPosPaths, extrPosDir, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_extr_pos.csv");
            std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> extrNegPaths;
            float_vec extrNegDir(-std::cos(EPSILON_OVER_TWO), -std::sin(EPSILON_OVER_TWO));
            infCyl2MC.traceRays(extrNegPaths, extrNegDir, infCyl2MC.M1.segments.back().second.y, h_inc, -infCyl2MC.M1.segments.back().second.y, "data_2mc_cyl_extr_neg.csv");
        }
    }



    // finite system
    if (FINITE_SYSTEM) {

        const TargetShape ts = TargetShape::CYLINDRICAL;
        const float_type f1(1), L(10), f2(1), da(0.0001), a_max(135 * DEG_TO_RAD);
        TwoMirrorConcentrator finCyl2MC;
        finCyl2MC = buildFin2MC(ts, f1, L, f2, da, a_max, "data_2mc_fin_cyl.csv");

        ParabolicMirror finCylPM;
        // gonna have to build a 2PM concentrator ...

    }





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

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////