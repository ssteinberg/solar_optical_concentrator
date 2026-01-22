#pragma once

#include "geometry/line_segment.h"

struct BoundingBox {
    float_vec min, max;

    BoundingBox();

    void fit(const LineSegment& l);

    bool intersect(const Ray& ray, HitInfo& minHit) const;
};
