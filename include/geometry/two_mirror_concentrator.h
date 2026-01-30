#pragma once

#include "geometry/geometry.h"
#include "geometry/mirror.h"
#include "geometry/polygon.h"

struct TwoMirrorConcentrator : Geometry {
    Mirror m1a, m1b, m2a, m2b, barrier;

    TwoMirrorConcentrator() : Geometry(Shape::CONSTRUCTED, Type::MIRROR) {}

    void reset();

    float_vec getCentre() const override;

    float_type getLength() const override;

    void buildFin(const bool& inv, const float_type& f1, const float_type& L, const float_type& f2, const float_type& da, const float_type& a_max, const float_type& w);

    void buildInf(const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max);


    bool calcCone(const Polygon& p, const float_vec& apex, HitInfo& h);

    void buildInfArb(const Polygon& p, const bool& inv, const float_type& L, const float_type& f, const float_vec& K_in, const float_type& dB, const float_type& B_max, const int& i, const int& j);


    void buildBarrier();

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;

    void writeTwoMirrorConcentrator(const std::string& filePath) const;

    [[nodiscard]] static float_type getAngularIntensityDistribution(float_type beta);
};
