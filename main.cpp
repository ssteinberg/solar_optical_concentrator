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

    // two mirrors
    Mirror M1, M2;

    // compute the line segments of the two mirrors
    void build(const TargetShape& targetShape, const bool& inverted, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max) {

        // incoming intensity
        const auto& I = [=](const float_type alpha) {
            if (targetShape == TargetShape::FLAT) return std::cos(alpha);
            else if (targetShape == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

        // outgoing intensity
        const auto& S = [=](const float_type beta) {
            if (targetShape == TargetShape::FLAT) return std::cos(beta);
            else if (targetShape == TargetShape::CYLINDRICAL) return float_type(1);
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
            const float_type dy = (inverted ? -1 : 1) * S(B) / I * dB;
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

// main
int main(int argc, char* argv[]) {



    // compute two-mirror concentrator coordinates
    std::cout << "Building two-mirror concentrator... ";
    const TargetShape targetShape = TargetShape::CYLINDRICAL;
    bool inverted;
    float_type L, f;
    const float_vec K_in(-1, 0);
    const float_type dB = 0.001;
    float_type B_max;
    if (targetShape == TargetShape::FLAT) {
        inverted = false;
        f = 0.2;
        L = 8 * f;
        B_max = 80 * DEG_TO_RAD;
    }
    else if (targetShape == TargetShape::CYLINDRICAL) {
        inverted = true;
        f = 0.5;
        L = 6 * f;
        B_max = 150 * DEG_TO_RAD;
    }
    TwoMirrorConcentrator TMC;
    TMC.build(targetShape, inverted, L, f, K_in, dB, B_max);
    std::cout << "Done." << std::endl;



    // ray trace
    std::cout << "Tracing rays... ";
    std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> TMCa, TMCb, TMCc;
    float_type extrRayDirX = -std::cos(EPSILON_OVER_TWO);
    float_type extrRayDirY = std::sin(EPSILON_OVER_TWO);
    float_vec posExtrRayDir(extrRayDirX, extrRayDirY);
    float_vec negExtrRayDir(extrRayDirX, -extrRayDirY);
    const float_type y_max = 5;
    float_type y = -y_max;
    const float_type increment = 0.01;
    while (y < y_max) {
        const Ray meanRay(float_vec(0, y), float_vec(-1, 0));
        TMC.traceRay(TMCa, meanRay);
        const Ray posExtrRay(float_vec(0, y), posExtrRayDir);
        TMC.traceRay(TMCb, posExtrRay);
        const Ray negExtrRay(float_vec(0, y), negExtrRayDir);
        TMC.traceRay(TMCc, negExtrRay);
        y += increment;
    }
    std::cout << "Done." << std::endl;



    // export data
    std::cout << "Writing to .csv... ";
    std::ofstream file;

    // two-mirror concentrator coordinates
    file.open("data_2mc.csv");
    int i = 0;
    int m = static_cast<int>(TMC.M1.segments.size());
    for (; i < m; ++i) {
        file << TMC.M1.segments[i].first.x << "," << TMC.M1.segments[i].first.y << ",";
        file << TMC.M2.segments[i].first.x << "," << TMC.M2.segments[i].first.y << "\n";
    }
    file << TMC.M1.segments[i - 1].second.x << "," << TMC.M1.segments[i - 1].second.y << ",";
    file << TMC.M2.segments[i - 1].second.x << "," << TMC.M2.segments[i - 1].second.y << "\n";
    file.close();

    // ray paths
    file.open("data_mean_paths.csv");
    for (auto rayPath : TMCa) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();
    file.open("data_pos_extr_paths.csv");
    for (auto rayPath : TMCb) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();
    file.open("data_neg_extr_paths.csv");
    for (auto rayPath : TMCc) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();

    // finished
    std::cout << "Done." << std::endl;



    // parabolic mirror
    ParabolicMirror PM;
    PM.build(f, 0.01, B_max);
    file.open("data_pm.csv");
    for (auto segment : PM.M.segments) file << segment.first.x << "," << segment.first.y << "\n";
    file.close();
    std::vector<std::tuple<float_vec, float_vec, float_vec>> PMa, PMb, PMc;
    float_type h = -1;
    while (h < 1) {
        const Ray meanRay(float_vec(0, h), float_vec(-1, 0));
        PM.traceRay(PMa, meanRay);
        const Ray posExtrRay(float_vec(0, h), posExtrRayDir);
        PM.traceRay(PMb, posExtrRay);
        const Ray negExtrRay(float_vec(0, h), negExtrRayDir);
        PM.traceRay(PMc, negExtrRay);
        h += 0.01;
    }
    file.open("data_pm_a.csv");
    for (auto path : PMa) {
        file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
        file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
        file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
    }
    file.close();
    file.open("data_pm_b.csv");
    for (auto path : PMb) {
        file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
        file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
        file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
    }
    file.close();
    file.open("data_pm_c.csv");
    for (auto path : PMc) {
        file << std::get<0>(path).x << "," << std::get<0>(path).y << ",";
        file << std::get<1>(path).x << "," << std::get<1>(path).y << ",";
        file << std::get<2>(path).x << "," << std::get<2>(path).y << "\n";
    }
    file.close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////