#pragma once

#include "types.h"
#include "ray.h"
#include "hitinfo.h"

struct Geometry {
    Shape shape;
    Type type;
    Geometry(const Shape& s, const Type& t) : shape(s), type(t) {}
    virtual float_vec getCentre() const = 0;
    virtual float_type getLength() const = 0;
    virtual bool intersect(const Ray& ray, HitInfo& hitInfo) const = 0;
    virtual Ray sampleMeanRay() const = 0;
    virtual std::pair<Ray, float_type> sampleDiffuseRay() const = 0;
    virtual std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeDiffuseRays() const = 0;
    virtual std::pair<std::vector<Ray>, std::vector<Ray>> generateExtremeInfiniteRays(const int& numRays) const = 0;
    virtual std::vector<Ray> generateFinalPlotRays() const = 0;
};
