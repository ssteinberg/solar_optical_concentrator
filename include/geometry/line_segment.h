#pragma once

#include "types.h"
#include "geometry/geometry.h"

struct LineSegment : Geometry {

    float_vec p1, p2, n1, n2;

    LineSegment(const Type &t, const float_vec &p1, const float_vec &p2, const float_vec &n1, const float_vec &n2): Geometry(Shape::FLAT, t), p1(p1), p2(p2), n1(n1), n2(n2) {}

    float_vec getCentre() const override;

    float_type getLength() const override;

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;

    void writeLineSegment(std::ofstream& file) const;
};