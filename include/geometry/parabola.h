#pragma once

#include "geometry/geometry.h"

struct Parabola : Geometry {
    float_type a, c;
    Parabola(const Type& t, const float_type& a, const float_type& c) : Geometry(Shape::PARABOLIC, t), a(a), c(c) {}
    float_vec getCentre() const override;

    float_type getLength() const override;

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;
};
