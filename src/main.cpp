////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// standard library
#include <filesystem>
#include <fstream>
#include <vector>

// project
#include "config.h"
#include "types.h"
#include "design.h"
#include "geometry/circle.h"
#include "geometry/ellipse.h"
#include "geometry/parabola.h"
#include "geometry/two_mirror_concentrator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void build_finite() {
#if !BUILD_FINITE
    return;
#endif

    // output
    std::string outputDataPath = "datafin/";
    if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
    std::ofstream file;

    // input
    const bool inv(false);
    const float_type f1(1), L(30), f2(1), da(0.000001), a_max(75 * DEG_TO_RAD), w(1), hl(f1 / 1000);
    file.open(outputDataPath + "input.csv");
    file << FLAT_SOURCE << "," << TARGET_SHAPE << "," << inv << "," << hl << "," << f1 << "," << L << "," << f2 << "," << da << "," << a_max * RAD_TO_DEG << "," << w << "\n";
    file.close();

    // two-mirror concentrator
    TwoMirrorConcentrator tmc;
    tmc.buildFin(inv, f1, L, f2, da, a_max, w);
    tmc.writeTwoMirrorConcentrator(outputDataPath);

    // design
    Design design;
    design.addGeometry(&tmc);
    LineSegment flatSource(Type::SOURCE, float_vec(-L, -hl), float_vec(-L, hl), float_vec(-1, 0), float_vec(-1, 0));
    Circle cylindricalSource(Type::SOURCE, float_vec(-L, 0), hl);
    LineSegment flatTarget(Type::TARGET, float_vec(0, -hl * w), float_vec(0, hl * w), float_vec(1, 0), float_vec(1, 0));
    Circle cylindricalTarget(Type::TARGET, float_vec(0, 0), hl);
    if (FLAT_SOURCE) {
        design.addGeometry(&flatSource);
        design.traceExtremeDiffuseRays(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) design.addGeometry(&flatTarget);
        else {
            cylindricalTarget.r /= PI;
            design.addGeometry(&cylindricalTarget);
        }
    }
    else {
        design.addGeometry(&cylindricalSource);
        design.traceExtremeDiffuseRays(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) {
            flatTarget.p1.y *= PI;
            flatTarget.p2.y *= PI;
            design.addGeometry(&flatTarget);
        }
        else design.addGeometry(&cylindricalTarget);
    }

    // final plot
    design.traceFinalPlotRays(outputDataPath);

    // phase space
    // design.tracePhaseSpace(outputDataPath, 100000);

    // hit report
    // const int numberOfDesigns = 10;
    // std::vector<Design> designs(numberOfDesigns);
    // std::vector<LineSegment> flatTargets;
    // std::vector<Circle> cylindricalTargets;
    // float_type increment = 2 * hl / numberOfDesigns;
    // if constexpr (FLAT_SOURCE && TARGET_SHAPE == Shape::CYLINDRICAL) increment /= PI;
    // else if constexpr (!FLAT_SOURCE && TARGET_SHAPE == Shape::FLAT) increment *= PI;
    // for (int i = 0; i < numberOfDesigns; ++i) {
    //     if constexpr (TARGET_SHAPE == Shape::FLAT) flatTargets.emplace_back(Type::TARGET, float_vec(0, -(i + 1) * increment), float_vec(0, (i + 1) * increment), float_vec(1, 0), float_vec(1, 0));
    //     else cylindricalTargets.emplace_back(Type::TARGET, float_vec(0, 0), (i + 1) * increment);
    // }
    // for (int i = 0; i < numberOfDesigns; ++i) {
    //     if (FLAT_SOURCE) designs[i].addGeometry(&flatSource);
    //     else designs[i].addGeometry(&cylindricalSource);
    //     if constexpr (TARGET_SHAPE == Shape::FLAT) designs[i].addGeometry(&flatTargets[i]);
    //     else designs[i].addGeometry(&cylindricalTargets[i]);
    //     designs[i].addGeometry(&tmc);
    // }
    // std::vector<float_vec> hitData;
    // const int numberOfTrials = 1;
    // const int numberOfRays = 10000;
    // for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitData(numberOfRays));
    // file.open(outputDataPath + "hitdata.csv");
    // for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
    // file.close();
}

void build_finite_ellipse() {
#if !BUILD_FINITE_ELLIPSE
    return;
#endif

    // output
    std::string outputDataPath = "datafinell/";
    if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
    std::ofstream file;

    // input
    const float_type f1(1), L(30), hl(f1 / 1000);
    file.open(outputDataPath + "input.csv");
    file << FLAT_SOURCE << "," << TARGET_SHAPE << "," << hl << "," << f1 << "," << L << "\n";
    file.close();

    // design
    Design design;
    Ellipse Elliot(Type::MIRROR, f1, L);
    design.addGeometry(&Elliot);
    LineSegment flatSource(Type::SOURCE, float_vec(-L, -hl), float_vec(-L, hl), float_vec(-1, 0), float_vec(-1, 0));
    Circle cylindricalSource(Type::SOURCE, float_vec(-L, 0), hl);
    LineSegment flatTarget(Type::TARGET, float_vec(0, -hl), float_vec(0, hl), float_vec(1, 0), float_vec(1, 0));
    Circle cylindricalTarget(Type::TARGET, float_vec(0, 0), hl);
    if (FLAT_SOURCE) {
        design.addGeometry(&flatSource);
        design.traceExtremeDiffuseRaysEllipse(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) design.addGeometry(&flatTarget);
        else {
            cylindricalTarget.r /= PI;
            design.addGeometry(&cylindricalTarget);
        }
    }
    else {
        design.addGeometry(&cylindricalSource);
        design.traceExtremeDiffuseRaysEllipse(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) {
            flatTarget.p1.y *= PI;
            flatTarget.p2.y *= PI;
            design.addGeometry(&flatTarget);
        }
        else design.addGeometry(&cylindricalTarget);
    }

    // phase space
    design.tracePhaseSpaceEllipse(outputDataPath, 100000);

    // hit report
    const int numberOfDesigns = 10;
    std::vector<Design> designs(numberOfDesigns);
    std::vector<LineSegment> flatTargets;
    std::vector<Circle> cylindricalTargets;
    float_type increment = 2 * hl / numberOfDesigns;
    if constexpr (FLAT_SOURCE && TARGET_SHAPE == Shape::CYLINDRICAL) increment /= PI;
    else if constexpr (!FLAT_SOURCE && TARGET_SHAPE == Shape::FLAT) increment *= PI;
    for (int i = 0; i < numberOfDesigns; ++i) {
        if constexpr (TARGET_SHAPE == Shape::FLAT) flatTargets.emplace_back(Type::TARGET, float_vec(0, -(i + 1) * increment), float_vec(0, (i + 1) * increment), float_vec(1, 0), float_vec(1, 0));
        else cylindricalTargets.emplace_back(Type::TARGET, float_vec(0, 0), (i + 1) * increment);
    }
    for (int i = 0; i < numberOfDesigns; ++i) {
        if (FLAT_SOURCE) designs[i].addGeometry(&flatSource);
        else designs[i].addGeometry(&cylindricalSource);
        if constexpr (TARGET_SHAPE == Shape::FLAT) designs[i].addGeometry(&flatTargets[i]);
        else designs[i].addGeometry(&cylindricalTargets[i]);
        designs[i].addGeometry(&Elliot);
    }
    std::vector<float_vec> hitData;
    const int numberOfTrials = 1;
    const int numberOfRays = 10000;
    for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitDataEllipse(numberOfRays));
    file.open(outputDataPath + "hitdata.csv");
    for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
    file.close();
}

void build_finite_parabola() {
#if !BUILD_FINITE_PARABOLA
    return;
#endif

    // output
    std::string outputDataPath = "datafinpara/";
    if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
    std::ofstream file;

    // input
    const float_type f1(1), L(4), f2(2), hl(f1 / 1000);
    file.open(outputDataPath + "input.csv");
    file << FLAT_SOURCE << "," << TARGET_SHAPE << "," << hl << "," << f1 << "," << L << "," << f2 << "\n";
    file.close();

    // design
    Design design;
    Parabola Param(Type::MIRROR1, 1 / (4 * f1), -f1 - L);
    Parabola Parker(Type::MIRROR2, -1 / (4 * f2), f2);
    design.addGeometry(&Param);
    design.addGeometry(&Parker);
    LineSegment flatSource(Type::SOURCE, float_vec(-L, -hl), float_vec(-L, hl), float_vec(-1, 0), float_vec(-1, 0));
    Circle cylindricalSource(Type::SOURCE, float_vec(-L, 0), hl);
    LineSegment flatTarget(Type::TARGET, float_vec(0, -hl), float_vec(0, hl), float_vec(1, 0), float_vec(1, 0));
    Circle cylindricalTarget(Type::TARGET, float_vec(0, 0), hl);
    if (FLAT_SOURCE) {
        design.addGeometry(&flatSource);
        design.traceExtremeDiffuseRaysParabola(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) design.addGeometry(&flatTarget);
        else {
            cylindricalTarget.r /= PI;
            design.addGeometry(&cylindricalTarget);
        }
    }
    else {
        design.addGeometry(&cylindricalSource);
        design.traceExtremeDiffuseRaysParabola(outputDataPath); // extreme rays
        if constexpr (TARGET_SHAPE == Shape::FLAT) {
            flatTarget.p1.y *= PI;
            flatTarget.p2.y *= PI;
            design.addGeometry(&flatTarget);
        }
        else design.addGeometry(&cylindricalTarget);
    }

    // phase space
    design.tracePhaseSpaceParabola(outputDataPath, 100000);

    // hit report
    // const int numberOfDesigns = 10;
    // std::vector<Design> designs(numberOfDesigns);
    // std::vector<LineSegment> flatTargets;
    // std::vector<Circle> cylindricalTargets;
    // float_type increment = 2 * hl / numberOfDesigns;
    // if constexpr (FLAT_SOURCE && TARGET_SHAPE == Shape::CYLINDRICAL) increment /= PI;
    // else if constexpr (!FLAT_SOURCE && TARGET_SHAPE == Shape::FLAT) increment *= PI;
    // for (int i = 0; i < numberOfDesigns; ++i) {
    //     if constexpr (TARGET_SHAPE == Shape::FLAT) flatTargets.emplace_back(Type::TARGET, float_vec(0, -(i + 1) * increment), float_vec(0, (i + 1) * increment), float_vec(1, 0), float_vec(1, 0));
    //     else cylindricalTargets.emplace_back(Type::TARGET, float_vec(0, 0), (i + 1) * increment);
    // }
    // for (int i = 0; i < numberOfDesigns; ++i) {
    //     if (FLAT_SOURCE) designs[i].addGeometry(&flatSource);
    //     else designs[i].addGeometry(&cylindricalSource);
    //     if constexpr (TARGET_SHAPE == Shape::FLAT) designs[i].addGeometry(&flatTargets[i]);
    //     else designs[i].addGeometry(&cylindricalTargets[i]);
    //     designs[i].addGeometry(&Param);
    //     designs[i].addGeometry(&Parker);
    // }
    // std::vector<float_vec> hitData;
    // const int numberOfTrials = 1;
    // const int numberOfRays = 10000;
    // for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceHitDataParabola(numberOfRays));
    // file.open(outputDataPath + "hitdata.csv");
    // for (const auto& p : hitData) file << p.x << "," << p.y << "\n";
    // file.close();
}

void build_infinite() {
#if !BUILD_INFINITE
    return;
#endif

    // output
    std::string outputDataPath = "datainf/";
    if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
    std::ofstream file;

    // input
    const bool inv(true);
    const float_type L(8), f(1), dB(0.00001), B_max(85 * DEG_TO_RAD), apertureSize(10);
    const float_vec K_in(-1, 0);
    file.open(outputDataPath + "input.csv");
    file << TARGET_SHAPE << "," << ELLIPTICAL_TARGET_X_RADIUS << "," << ELLIPTICAL_TARGET_Y_RADIUS << "," << inv << "," << L << "," << f << "," << dB << "," << B_max * RAD_TO_DEG << "," << apertureSize << "\n";
    file.close();

    // two-mirror concentrator
    TwoMirrorConcentrator tmc;
    tmc.buildInf(inv, L, f, K_in, dB, B_max, apertureSize);
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

void build_infinite_arbitrary() {
#if !BUILD_INFINITE_ARBITRARY
    return;
#endif

    // output
    std::string outputDataPath = "datainfarb/";
    if (!std::filesystem::exists(outputDataPath)) std::filesystem::create_directory(outputDataPath);
    std::ofstream file;

    // input
    const bool inv(true);
    const float_type L(3), f(0.5), dB(0.00001), B_max(85 * DEG_TO_RAD);
    const float_vec K_in(-1, 0);
    file.open(outputDataPath + "input.csv");
    file << EPSILON << "," << inv << "," << L << "," << f << "," << dB << "," << B_max * RAD_TO_DEG << "\n";
    file.close();

    // infinite source
    LineSegment source(Type::SOURCE, float_vec(0, -2), float_vec(0, 2), K_in, K_in);

    // target
    const std::vector<float_vec> targetVertices = {float_vec(0, -EPSILON_OVER_TWO), float_vec(0, EPSILON_OVER_TWO)}; // ideal
    // const std::vector<float_vec> targetVertices = {float_vec(1, 1), float_vec(-2, 0), float_vec(0, -1)}; // triangle
    const Polygon target(Type::TARGET, targetVertices);

    // two-mirror concentrator
    TwoMirrorConcentrator tmc;
    tmc.buildInfArb(target, inv, L, f, K_in, dB, B_max, 6, 2);
    tmc.writeTwoMirrorConcentrator(outputDataPath);

    // design
    Design d;
    d.addGeometry(&tmc);
    d.addGeometry(&source);

    // trace extreme rays
    d.traceExtremeInfiniteRays(outputDataPath, 50);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
        build_finite();
        build_finite_ellipse();
        build_finite_parabola();
        build_infinite();
        build_infinite_arbitrary();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

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
if constexpr (TARGET_SHAPE == Shape::FLAT) design.addGeometry(&flatTarget);
else if (CYLINDRICAL_TARGET) design.addGeometry(&cylindricalTarget);
auto phaseSpaceData = design.tracePhaseSpace(10000);
file.open(outputDataPath + "phase.csv");
for (auto data : phaseSpaceData) file << data[0] << "," << data[1] << "," << data[2] << "\n";
file.close();

*/

/*

// flat target
const LineSegment target(Type::TARGET, float_vec(0, -EPSILON_OVER_TWO), float_vec(0, EPSILON_OVER_TWO), float_vec(1, 0), float_vec(1, 0));

// flat target
const float_vec tp1(0, -EPSILON_OVER_TWO), tp2(0, EPSILON_OVER_TWO);
const float_vec tn(1, 0);
Mirror t; t.addSegment(tp1, tp2, tn, tn);
const std::vector<float_vec> tpts = {tp1, tp2};

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

/*
std::vector<std::pair<float_type, std::vector<float_type>>> hitData;
const int numberOfTrials = 3;
const int numberOfRays = 10000;
for (const auto& d : designs) for (int i = 0; i < numberOfTrials; ++i) hitData.emplace_back(d.target->getLength() / d.source->getLength(), d.traceDetailedHitData(numberOfRays));
file.open(outputDataPath + "hitdata.csv");
for (const auto& p : hitData) file << p.first << "," << p.second[0] << "," << p.second[1] << "," << p.second[2] << "," << p.second[3] << "\n";
file.close();
*/
