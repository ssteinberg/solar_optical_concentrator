////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//  compile: "g++-15 main.cpp"
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

// global constants
constexpr float PI = 3.14159265358979f;
constexpr float DEG_TO_RAD = PI / 180.0f;
constexpr float RAD_TO_DEG = 180.0f / PI;
float EPSILON = 0.01f;
float EPSILON_OVER_TWO = EPSILON / 2.0f;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// ray
struct Ray {
    float2 o, d;
    Ray(const float2& o, const float2& d) : o(o), d(d) {}
};

// hit info
struct HitInfo {
    float2 P, N;
};

// target shape
enum struct TargetShape {
    FLAT,
    CYLINDRICAL
};

// mirror
struct Mirror {

    // coordinates and normals of the line segments of the mirror
    std::vector<std::pair<float2, float2>> segments;
    std::vector<std::pair<float2, float2>> normals;

    // add coordinates and normals of a line segment to the mirror
    void addSegment(const float2& p1, const float2& p2, const float2& n1, const float2& n2) {
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
            float2 segmentDirection = segments[i].second - segments[i].first;
            float denom = cross(segmentDirection, ray.d);
            if (denom == 0) continue;

            // intersection
            float t = cross((ray.o - segments[i].first), ray.d) / denom;
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
    void build(const TargetShape& targetShape, const bool& inverted, const float& L, const float& f, const float2& K_in, const float& dB, const float& B_max) {

        // outgoing intensity to target
        const auto& S = [=](const float beta) {
            if (targetShape == TargetShape::FLAT) return std::cos(beta);
            else if (targetShape == TargetShape::CYLINDRICAL) return 1.0f;
            else return 0.0f;
        };

        // initial conditions
        float2 p1(-L, 0.0f);
        float2 p2(f, 0.0f);
        float2 K_int = p2 - p1;
        float R = length(K_int);
        K_int /= R;
        float2 K_out = -p2;
        float r = length(K_out);
        K_out /= r;
        float2 n1 = K_int - K_in;
        float2 n2 = K_out - K_int;

        // numerical integration
        float B = 0.0f;
        while (B < B_max) {

            // step in beta
            float B_new = B + dB;

            // precomputation
            float sinB = std::sin(B);
            float cosB = std::cos(B);
            float d = R - 2.0f * (L + f);

            // (7)
            float dy = (inverted ? -1.0f : 1.0f) * S(B) * dB;
            float dx = dy * (p1.y - r * sinB) / (r * (cosB - 1.0f) + 2.0f * (L + f));
            float2 p1_new(p1.x + dx, p1.y + dy);

            // (11)
            float dr = r * ((r + d) * sinB - p1.y * cosB) / ((r + d) * cosB + p1.y * sinB - r - R) * dB;
            float r_new = r + dr;
            float2 p2_new = r_new * float2(std::cos(B_new), std::sin(B_new));

            // new directions
            float2 K_int_new = p2_new - p1_new;
            float R_new = length(K_int_new);
            K_int_new /= R_new;
            float2 K_out_new = -p2_new / r_new;
            float2 n1_new = K_int_new - K_in;
            float2 n2_new = K_out_new - K_int_new;

            // store coordinates and normals
            M1.addSegment(p1, p1_new, n1, n1_new);
            M1.addSegment(float2(p1.x, -p1.y), float2(p1_new.x, -p1_new.y), float2(n1.x, -n1.y), float2(n1_new.x, -n1_new.y));
            M2.addSegment(p2, p2_new, n2, n2_new);
            M2.addSegment(float2(p2.x, -p2.y), float2(p2_new.x, -p2_new.y), float2(n2.x, -n2.y), float2(n2_new.x, -n2_new.y));

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
    void traceRay(std::vector<std::tuple<float2, float2, float2, float2>>& rayPaths, const Ray& ray1) {

        // check intersection with mirror 1
        if (HitInfo hit1; M1.intersect(hit1, ray1)) {

            // check intersection with mirror 2
            const Ray ray2(hit1.P, normalize(ray1.d - 2 * dot(hit1.N, ray1.d) * hit1.N));
            if (HitInfo hit2; M2.intersect(hit2, ray2)) {

                // compute point through target and store path of ray
                rayPaths.emplace_back(ray1.o, hit1.P, hit2.P, hit2.P + 10.0f * normalize(ray2.d - 2 * dot(hit2.N, ray2.d) * hit2.N));
            }
        }
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// main
int main(int argc, char* argv[]) {



    // compute two-mirror concentrator coordinates
    const TargetShape targetShape = TargetShape::CYLINDRICAL;
    bool inverted;
    float L;
    const float f = 1.0f;
    const float2 K_in(-1.0f, 0.0f);
    float dB = 0.001f;
    float B_max;
    if (targetShape == TargetShape::FLAT) {
        inverted = false;
        L = 8.0f * f;
        B_max = 80.0f * DEG_TO_RAD;
    }
    else if (targetShape == TargetShape::CYLINDRICAL) {
        inverted = true;
        L = 6.0f * f;
        B_max = 150.0f * DEG_TO_RAD;
    }
    TwoMirrorConcentrator TMC;
    TMC.build(targetShape, inverted, L, f, K_in, dB, B_max);



    // ray trace
    std::vector<std::tuple<float2, float2, float2, float2>> rayPaths;

    // extreme rays
    float extremeRayDirectionX = -std::cos(EPSILON_OVER_TWO);
    float extremeRayDirectionY = std::sin(EPSILON_OVER_TWO);
    float2 extremeRayDirectionPositive(extremeRayDirectionX, extremeRayDirectionY);
    float2 extremeRayDirectionNegative(extremeRayDirectionX, -extremeRayDirectionY);

    // iterate over y
    const float y_max = 3.0f;
    float y = -y_max;
    const float increment = 0.1f;
    while (y < y_max) {

        // trace positive extreme ray
        const Ray extremeRayPositive(float2(0.0f, y), extremeRayDirectionPositive);
        TMC.traceRay(rayPaths, extremeRayPositive);

        // trace negative extreme ray
        const Ray extremeRayNegative(float2(0.0f, y), extremeRayDirectionNegative);
        TMC.traceRay(rayPaths, extremeRayNegative);

        // update y
        y += increment;
    }



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

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////