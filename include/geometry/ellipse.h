#pragma once

#include "geometry/geometry.h"

struct Ellipse : Geometry {
    float_type f1, L, saa, sab;
    Ellipse(const Type& t, const float_type& f1, const float_type& L) : Geometry(Shape::ELLIPTICAL, t), f1(f1), L(L), saa(f1 + L / 2), sab(std::sqrt(f1 * (f1 + L))) {}
    float_vec getCentre() const override;

    float_type getLength() const override;

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;
};
