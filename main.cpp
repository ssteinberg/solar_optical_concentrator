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
    void build1(const TargetShape& targetShape, const bool& inverted, const float_type& L, const float_type& f, const float_type& B_max) {

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

        // parameters
        float_vec K_in(-1, 0);
        float_type dB = 0.01;

        // initial conditions
        float_vec p(-L, 0), pp(f, 0);
        float_vec n(1, 0), np(-1, 0);
        float_type R(L + f), r(f);

        // numerical integration
        float_type B(0);
        while (std::abs(B) < B_max) {

            // step in beta
            float_type B_new = B + dB;

            // precomputation
            float_type sinB = std::sin(B);
            float_type cosB = std::cos(B);
            float_type d = R - 2 * (L + f);

            // (7)
            float_type dy = (inverted ? -1 : 1) * S(B) * dB;
            float_type dx = dy * (p.y - r * sinB) / (r * (cosB - 1) + 2 * (L + f));
            float_vec p_new(p.x + dx, p.y + dy);

            // (11)
            float_type dr = r * ((r + d) * sinB - p.y * cosB) / ((r + d) * cosB + p.y * sinB - r - R) * dB;
            float_type r_new = r + dr;
            float_vec pp_new = r_new * float_vec(std::cos(B_new), std::sin(B_new));

            // new directions
            float_vec K_int = pp_new - p_new;
            float_type R_new = length(K_int);
            K_int /= R_new;
            float_vec K_out = -pp_new / r_new;
            float_vec n_new = normalize(K_int - K_in);
            float_vec np_new = normalize(K_out - K_int);

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

    // compute the line segments of the two mirrors
    void build2(const TargetShape& targetShape, const bool& inverted, const float_type& L, const float_type& f, const float_type& B_max) {

        // incoming intensity to M1
        const auto& I = [=](const float_type alpha) {
            if (targetShape == TargetShape::FLAT) return std::cos(alpha);
            else if (targetShape == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

        // outgoing intensity to target
        const auto& S = [=](const float_type beta) {
            if (targetShape == TargetShape::FLAT) return std::cos(beta);
            else if (targetShape == TargetShape::CYLINDRICAL) return float_type(1);
            else return float_type(0);
        };

        // parameters
        const float_vec K_in(-std::cos(EPSILON), -std::sin(EPSILON));
        const float_type dy = 0.01 * std::cos(EPSILON_OVER_TWO);

        // initial conditions
        float_type x(-L), y(0), xp(f), yp(0), R(L + f), r(f);
        float_vec n(std::cos(EPSILON_OVER_TWO), std::sin(EPSILON_OVER_TWO)), np(-1, 0);
        float_type dx_dy(-n.y / n.x);

        // numerical integration
        float_type B(0), Br(0);
        while (std::abs(B) < B_max) {

            const float_type y_new = y + dy;
            const float_type dx = dx_dy * (y_new - y);
            const float_type x_new = x + dx;

            const float_type Br_new = cross(K_in, float_vec(x_new + L, y_new));
            const float_type dB = (inverted ? 1 : -1) * I(Br) / S(B) * (Br_new - Br); // * sign(dy)
            const float_type B_new = B + dB;

            const float_vec step(np.y, -np.x);
            const float_vec p_hat(std::cos(B), std::sin(B));
            const float_type sina = cross(p_hat, step);

            const float_vec p2_new = float_vec(xp, yp) + r * std::sin(dB) / sina * step;
            const float_type xp_new = p2_new.x;
            const float_type yp_new = p2_new.y;

            const float_type r_new = sqrt(xp_new * xp_new + yp_new * yp_new);
            const float_type R_new = sqrt((xp_new - x_new) * (xp_new - x_new) + (yp_new - y_new) * (yp_new - y_new));
            const float_vec K_int = float_vec(xp_new - x_new, yp_new - y_new) / R_new;
            const float_vec K_out = -float_vec(std::cos(B_new), std::sin(B_new));

            float_vec n_new = K_int - K_in;
            float_vec np_new = K_out - K_int;
            const float_type dx_dy_new = -n_new.y / n_new.x;
            n_new = normalize(n_new);
            np_new = normalize(np_new);

            M1.addSegment(float_vec(x, y), float_vec(x_new, y_new), n, n_new);
            M2.addSegment(float_vec(xp, yp), p2_new, np, np_new);

            x = x_new;
            y = y_new;
            xp = xp_new;
            yp = yp_new;
            R = R_new;
            r = r_new;
            dx_dy = dx_dy_new;
            n = n_new;
            np = np_new;
            B = B_new;
            Br = Br_new;
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

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// main
int main(int argc, char* argv[]) {



    // compute two-mirror concentrator coordinates
    std::cout << "Building two-mirror concentrator... ";
    const TargetShape targetShape = TargetShape::CYLINDRICAL;
    bool inverted;
    float_type L;
    const float_type f = 1;
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
    TMC.build1(targetShape, inverted, L, f, B_max);
    // TMC.build2(targetShape, inverted, L, f, B_max);
    std::cout << "Done." << std::endl;



    // ray trace
    std::cout << "Tracing rays... ";
    std::vector<std::tuple<float_vec, float_vec, float_vec, float_vec>> meanPaths, posExtrPaths, negExtrPaths;
    float_type extrRayDirX = -std::cos(EPSILON_OVER_TWO);
    float_type extrRayDirY = std::sin(EPSILON_OVER_TWO);
    float_vec posExtrRayDir(extrRayDirX, extrRayDirY);
    float_vec negExtrRayDir(extrRayDirX, -extrRayDirY);
    const float_type y_max = 5;
    float_type y = -y_max;
    const float_type increment = 0.01;
    while (y < y_max) {
        const Ray meanRay(float_vec(0, y), float_vec(-1, 0));
        TMC.traceRay(meanPaths, meanRay);
        const Ray posExtrRay(float_vec(0, y), posExtrRayDir);
        TMC.traceRay(posExtrPaths, posExtrRay);
        const Ray negExtrRay(float_vec(0, y), negExtrRayDir);
        TMC.traceRay(negExtrPaths, negExtrRay);
        y += increment;
    }
    std::cout << "Done." << std::endl;



    // export data
    std::cout << "Writing to .csv... ";
    std::ofstream file;

    // two-mirror concentrator coordinates
    file.open("data_coords.csv");
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
    for (auto rayPath : meanPaths) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();
    file.open("data_pos_extr_paths.csv");
    for (auto rayPath : posExtrPaths) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();
    file.open("data_neg_extr_paths.csv");
    for (auto rayPath : negExtrPaths) {
        file << std::get<0>(rayPath).x << "," << std::get<0>(rayPath).y << ",";
        file << std::get<1>(rayPath).x << "," << std::get<1>(rayPath).y << ",";
        file << std::get<2>(rayPath).x << "," << std::get<2>(rayPath).y << ",";
        file << std::get<3>(rayPath).x << "," << std::get<3>(rayPath).y << "\n";
    }
    file.close();

    // finished
    std::cout << "Done." << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////