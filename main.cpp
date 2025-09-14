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

// linear algebra
#include "linalg.h"
using namespace linalg::aliases;

// floating-point precision
using float_type = double;
constexpr int NUMBER_OF_DIMENSIONS = 2;
typedef linalg::vec<float_type, NUMBER_OF_DIMENSIONS> float_vec;

// global constants
constexpr float_type PI = 3.14159265358979;
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

        // outgoing intensity to target
        const auto& S = [=](const float_type beta) {
            if (targetShape == TargetShape::FLAT) return std::cos(beta);
            else if (targetShape == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

        // initial conditions
        float_vec p1(-L, 0);
        float_vec p2(f, 0);
        float_vec K_int = p2 - p1;
        float_type R = length(K_int);
        K_int /= R;
        float_vec K_out = -p2;
        float_type r = length(K_out);
        K_out /= r;
        float_vec n1 = K_int - K_in;
        float_vec n2 = K_out - K_int;

        // numerical integration
        float_type B = 0;
        while (B < B_max) {

            // step in beta
            float_type B_new = B + dB;

            // precomputation
            float_type sinB = std::sin(B);
            float_type cosB = std::cos(B);
            float_type d = R - 2 * (L + f);

            // (7)
            float_type dy = (inverted ? -1 : 1) * S(B) * dB;
            float_type dx = dy * (p1.y - r * sinB) / (r * (cosB - 1) + 2 * (L + f));
            float_vec p1_new(p1.x + dx, p1.y + dy);

            // (11)
            float_type dr = r * ((r + d) * sinB - p1.y * cosB) / ((r + d) * cosB + p1.y * sinB - r - R) * dB;
            float_type r_new = r + dr;
            float_vec p2_new = r_new * float_vec(std::cos(B_new), std::sin(B_new));

            // new directions
            float_vec K_int_new = p2_new - p1_new;
            float_type R_new = length(K_int_new);
            K_int_new /= R_new;
            float_vec K_out_new = -p2_new / r_new;
            float_vec n1_new = K_int_new - K_in;
            float_vec n2_new = K_out_new - K_int_new;

            // store coordinates and normals
            M1.addSegment(p1, p1_new, n1, n1_new);
            M1.addSegment(float_vec(p1.x, -p1.y), float_vec(p1_new.x, -p1_new.y), float_vec(n1.x, -n1.y), float_vec(n1_new.x, -n1_new.y));
            M2.addSegment(p2, p2_new, n2, n2_new);
            M2.addSegment(float_vec(p2.x, -p2.y), float_vec(p2_new.x, -p2_new.y), float_vec(n2.x, -n2.y), float_vec(n2_new.x, -n2_new.y));

            // update
            p1 = p1_new;
            p2 = p2_new;
            K_int = K_int_new;
            K_out = K_out_new;
            n1 = n1_new;
            n2 = n2_new;
            r = r_new;
            R = R_new;
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
                const float_type len = 10;
                rayPaths.emplace_back(ray1.o, hit1.P, hit2.P, hit2.P + len * normalize(ray2.d - 2 * dot(hit2.N, ray2.d) * hit2.N));
            }
        }
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// main
int main(int argc, char* argv[]) {

    std::cout << "Building two-mirror concentrator ... ";

    // compute two-mirror concentrator coordinates
    const TargetShape targetShape = TargetShape::CYLINDRICAL;
    bool inverted;
    float_type L;
    const float_type f = 1;
    const float_vec K_in(-1, 0);
    float_type dB = 0.001;
    float_type B_max;
    if (targetShape == TargetShape::FLAT) {
        inverted = false;
        L = 8 * f;
        B_max = 80 * DEG_TO_RAD;
    }
    else if (targetShape == TargetShape::CYLINDRICAL) {
        inverted = true;
        L = 6 * f;
        B_max = 150 * DEG_TO_RAD;
    }
    TwoMirrorConcentrator TMC;
    TMC.build(targetShape, inverted, L, f, K_in, dB, B_max);

    std::cout << "Done." << std::endl;
    std::cout << "Tracing rays ... ";

    // ray trace
    std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> rayPaths;

    // extreme rays
    float_type extremeRayDirectionX = -std::cos(EPSILON_OVER_TWO);
    float_type extremeRayDirectionY = std::sin(EPSILON_OVER_TWO);
    float_vec extremeRayDirectionPositive(extremeRayDirectionX, extremeRayDirectionY);
    float_vec extremeRayDirectionNegative(extremeRayDirectionX, -extremeRayDirectionY);

    // iterate over y
    const float_type y_max = 3;
    float_type y = -y_max;
    const float_type increment = 0.1;
    while (y < y_max) {

        // trace positive extreme ray
        const Ray extremeRayPositive(float_vec(0, y), extremeRayDirectionPositive);
        TMC.traceRay(rayPaths, extremeRayPositive);

        // trace negative extreme ray
        const Ray extremeRayNegative(float_vec(0, y), extremeRayDirectionNegative);
        TMC.traceRay(rayPaths, extremeRayNegative);

        // update y
        y += increment;
    }

    std::cout << "Done." << std::endl;
    std::cout << "Writing to .csv ... ";

    // export data
    std::ofstream file;

    // two-mirror concentrator coordinates
    file.open("data_coords.csv");
    int i = 0;
    int m = static_cast<int>(TMC.M1.segments.size());
    for (; i < m; ++i) {
        file << TMC.M1.segments[i].first.x << "," << TMC.M1.segments[i].first.y << ",";
        file << TMC.M1.segments[i].second.x << "," << TMC.M1.segments[i].second.y << ",";
        file << TMC.M2.segments[i].first.x << "," << TMC.M2.segments[i].first.y << ",";
        file << TMC.M2.segments[i].second.x << "," << TMC.M2.segments[i].second.y << "\n";
    }
    file.close();

    // ray paths
    file.open("data_rays.csv");
    for (auto rayPath : rayPaths) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();

    std::cout << "Done." << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////