#pragma once

#include "geometry/geometry.h"
#include "geometry/line_segment.h"

struct Polygon : Geometry {
    std::vector<float_vec> vertices;
    std::vector<LineSegment> segments;

    Polygon(const Type& type, const std::vector<float_vec>& vertices);

    float_vec getCentre() const override;

    float_type getLength() const override;

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;

    void writePolygon(std::ofstream& file) const;
};
