#pragma once

#include "geometry/geometry.h"
#include "geometry/line_segment.h"
#include "tree.h"

struct Mirror : Geometry {
    std::vector<LineSegment> segments;
    Tree tree;

    Mirror() : Geometry(Shape::CONSTRUCTED, Type::MIRROR) {}

    void reset();

    float_vec getCentre() const override;

    float_type getLength() const override;

    void addSegment(const float_vec& p1, const float_vec& p2, const float_vec& n1, const float_vec& n2);

    void buildTree();

    bool intersect(const Ray& ray, HitInfo& hitInfo) const override;

    Ray sampleMeanRay() const override;

    std::pair<Ray, float_type> sampleDiffuseRay() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const override;

    std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const override;

    std::vector<Ray> generateFinalPlotRays() const override;

    void writeMirror(std::ofstream& file) const;
};
