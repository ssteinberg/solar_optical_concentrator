#pragma once

#include "geometry/geometry.h"

struct Circle : Geometry {
    float_vec c;
    float_type r;

    Circle(const Type &t, const float_vec &c, const float_type &r) : Geometry(Shape::CYLINDRICAL, t), c(c), r(r) {}

    float_vec getCentre() const override;

    float_type getLength() const override;

    bool intersect(const Ray &ray, HitInfo &hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray> > generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray> > generateExtremeInfiniteRays(const int &numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;
};
