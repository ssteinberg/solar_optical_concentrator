#include "geometry/mirror.h"
#include "pcg32.h"

void Mirror::reset() { segments.clear(); }

float_vec Mirror::getCentre() const { return float_vec(0, 0); }

float_type Mirror::getLength() const {
    float_type l(0);
    for (const auto& s : segments) l += s.getLength();
    return l;
}

void Mirror::addSegment(const float_vec &p1, const float_vec &p2, const float_vec &n1, const float_vec &n2) {
    segments.emplace_back(Type::MIRROR, p1, p2, n1, n2);
}

void Mirror::buildTree() {
    tree = Tree(segments);
}

bool Mirror::intersect(const Ray &ray, HitInfo &hitInfo) const {
    float_type length = 0;
    for (auto s : segments) {
        length += s.getLength();
        if (s.intersect(ray, hitInfo)) {
            hitInfo.ml = length;
            return true;
        }
    }
    return false;
    // return tree.intersect(ray, hitInfo);
}

Ray Mirror::sampleMeanRay() const {
    return segments[PCG32::rand() * segments.size()].sampleMeanRay();
}

std::pair<Ray, float_type> Mirror::sampleDiffuseRay() const {
    return segments[PCG32::rand() * segments.size()].sampleDiffuseRay();
}

std::pair<std::vector<Ray>, std::vector<Ray>> Mirror::generateExtremeDiffuseRays() const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::pair<std::vector<Ray>, std::vector<Ray>> Mirror::generateExtremeInfiniteRays(const int &numRays) const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::vector<Ray> Mirror::generateFinalPlotRays() const {
    return std::vector<Ray>();
}

void Mirror::writeMirror(std::ofstream &file) const {
    for (auto s : segments) {
        s.writeLineSegment(file);
    }
}
