#include "geometry/line_segment.h"
#include "pcg32.h"

#include <fstream>


float_vec LineSegment::getCentre() const { return (p1 + p2) / 2; }

float_type LineSegment::getLength() const { return length(p1 - p2); }

bool LineSegment::intersect(const Ray &ray, HitInfo &hitInfo) const {
    if (IGNORE_SOURCE && type == Type::SOURCE) return false;
    else if (dot(n1, -ray.d) <= 0 && type == Type::TARGET ) return false;
    else {
        const float_vec segDir = p2 - p1;
        const float_type rayDir_x_segDir = cross(ray.d, segDir);
        if (rayDir_x_segDir == 0) return false;
        const float_type s = cross((p1 - ray.o), ray.d) / rayDir_x_segDir;
        const float_type t = cross((p1 - ray.o), segDir) / rayDir_x_segDir;
        if (0 <= s && s <= 1 && 0 < t) {
            hitInfo.l = t;
            hitInfo.p = ray.o + hitInfo.l * ray.d;
            hitInfo.n = normalize((1 - s) * n1 + s * n2);
            hitInfo.t = type;
            return true;
        }
        return false;
    }
}

Ray LineSegment::sampleMeanRay() const {
    const float_type s = PCG32::rand();
    const float_vec p = (1 - s) * p1 + s * p2;
    const float_vec n = (1 - s) * n1 + s * n2;
    return Ray(p + (DOINKING ? DOINK * n : float_vec(0, 0)), n);
}

std::pair<Ray, float_type> LineSegment::sampleDiffuseRay() const {
    const Ray meanRay = sampleMeanRay();
    const float_vec n = meanRay.d;
    const float_type theta = std::asin(2 * PCG32::rand() - 1);
    const float_type cos0 = std::cos(theta);
    const float_type sin0 = std::sin(theta);
    const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
    return std::pair(Ray(meanRay.o, rotDir), theta * RAD_TO_DEG);
}

std::pair<std::vector<Ray>, std::vector<Ray>> LineSegment::generateExtremeDiffuseRays() const {
    std::vector<Ray> p1ExtrRays, p2ExtrRays;
    for (int degrees = -90; degrees <= 90; degrees += 2) {
        const float_type theta = degrees * DEG_TO_RAD;
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec p1RotDir(n1.x * cos0 - n1.y * sin0, n1.x * sin0 + n1.y * cos0);
        const float_vec p2RotDir(n2.x * cos0 - n2.y * sin0, n2.x * sin0 + n2.y * cos0);
        p1ExtrRays.emplace_back(p1 + (DOINKING ? DOINK * n1 : float_vec(0, 0)), p1RotDir);
        p2ExtrRays.emplace_back(p2 + (DOINKING ? DOINK * n2 : float_vec(0, 0)), p2RotDir);
    }
    return std::pair(p1ExtrRays, p2ExtrRays);
}

std::pair<std::vector<Ray>, std::vector<Ray>> LineSegment::generateExtremeInfiniteRays(const int &numRays) const {
    std::vector<Ray> posExtrRays, negExtrRays;
    const float_type inc(1 / static_cast<float_type>(numRays));
    float_type l(0);
    while (l <= 1) {
        const float_vec p((1 - l) * p1 + l * p2);
        const float_vec n((1 - l) * n1 + l * n2);
        const float_type cosPos(std::cos(EPSILON_OVER_TWO)), sinPos(std::sin(EPSILON_OVER_TWO));
        const float_type cosNeg(std::cos(-EPSILON_OVER_TWO)), sinNeg(std::sin(-EPSILON_OVER_TWO));
        const float_vec posRotDir(n.x * cosPos - n.y * sinPos, n.x * sinPos + n.y * cosPos);
        const float_vec negRotDir(n.x * cosNeg - n.y * sinNeg, n.x * sinNeg + n.y * cosNeg);
        posExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), posRotDir);
        negExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), negRotDir);
        l += inc;
    }
    return std::pair(posExtrRays, negExtrRays);
}

std::vector<Ray> LineSegment::generateFinalPlotRays() const {
    std::vector<Ray> rays;
    for (int degrees = -90; degrees <= 90; degrees += 15) {
        const float_type theta = degrees * DEG_TO_RAD;
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec p1RotDir(n1.x * cos0 - n1.y * sin0, n1.x * sin0 + n1.y * cos0);
        rays.emplace_back(getCentre() + (DOINKING ? DOINK * n1 : float_vec(0, 0)), p1RotDir);
    }
    return rays;
}

void LineSegment::writeLineSegment(std::ofstream &file) const {
    file << p1.x << "," << p1.y << "," << p2.x << "," << p2.y << "\n";
}
