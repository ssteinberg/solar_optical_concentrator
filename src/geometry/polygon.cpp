#include "geometry/polygon.h"
#include "pcg32.h"

#include <iostream>

Polygon::Polygon(const Type &type, const std::vector<float_vec> &vertices): Geometry(Shape::CONSTRUCTED, type), vertices(vertices) {
    const int numVertices(static_cast<int>(vertices.size()));
    if (numVertices < 2) std::cout << "Not a polygon." << std::endl;
    for (int i = 0; i < numVertices - 1; ++i) {
        const float_vec p1(vertices[i]), p2(vertices[i + 1]);
        const float_vec l(normalize(p2 - p1)), n(l.y, -l.x);
        segments.emplace_back(type, p1, p2, n, n);
    }
    if (numVertices > 2) {
        const float_vec p1(vertices.back()), p2(vertices.front());
        const float_vec l(normalize(p2 - p1)), n(l.y, -l.x);
        segments.emplace_back(type, p1, p2, n, n);
    }
}

float_vec Polygon::getCentre() const {
    float_vec sum(0, 0);
    for (const auto& s : segments) sum += s.p1;
    return sum / static_cast<int>(segments.size());
}

float_type Polygon::getLength() const {
    float_type l(0);
    for (const auto& s : segments) l += s.getLength();
    return l;
}

bool Polygon::intersect(const Ray &ray, HitInfo &hitInfo) const {
    for (auto s : segments) if (s.intersect(ray, hitInfo)) return true;
    return false;
}

Ray Polygon::sampleMeanRay() const {
    return segments[PCG32::rand() * segments.size()].sampleMeanRay();
}

std::pair<Ray, float_type> Polygon::sampleDiffuseRay() const {
    return segments[PCG32::rand() * segments.size()].sampleDiffuseRay();
}

std::pair<std::vector<Ray>, std::vector<Ray>> Polygon::generateExtremeDiffuseRays() const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::pair<std::vector<Ray>, std::vector<Ray>> Polygon::generateExtremeInfiniteRays(const int &numRays) const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::vector<Ray> Polygon::generateFinalPlotRays() const {
    return std::vector<Ray>();
}

void Polygon::writePolygon(std::ofstream &file) const {
    for (auto s : segments) s.writeLineSegment(file);
}

