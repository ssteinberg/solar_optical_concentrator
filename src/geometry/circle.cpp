#include "geometry/circle.h"
#include "pcg32.h"

float_vec Circle::getCentre() const { return c; }

float_type Circle::getLength() const { return PI * 2 * r; }

bool Circle::intersect(const Ray &ray, HitInfo &hitInfo) const {
    if (IGNORE_SOURCE && type == Type::SOURCE) return false;
    else {
        const float_vec oc = c - ray.o;
        const float_type hyp2 = dot(oc, oc);
        if (hyp2 < r * r) return false;
        const float_type hyp = std::sqrt(hyp2);
        const float_vec ocDir = oc / hyp;
        const float_type cos0 = dot(ray.d, ocDir);
        if (cos0 <= 0) return false;
        const float_type adj = hyp * cos0;
        const float_type d2 = hyp2 - adj * adj;
        if (d2 > r * r) return false;
        if (d2 == r * r) {
            hitInfo.l = adj;
            hitInfo.p = ray.o + hitInfo.l * ray.d;
        }
        else if (d2 < r * r) {
            hitInfo.l = adj - std::sqrt(r * r - d2);
            hitInfo.p = ray.o + hitInfo.l * ray.d;
        }
        hitInfo.n = normalize(hitInfo.p - c);
        hitInfo.t = type;
        return true;
    }
}

Ray Circle::sampleMeanRay() const {
    const float_type theta = PCG32::rand() * 2 * PI;
    const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
    const float_vec n = (p - c) / r;
    return Ray(p + (DOINKING ? DOINK * n : float_vec(0, 0)), n);
}

std::pair<Ray, float_type> Circle::sampleDiffuseRay() const {
    const Ray meanRay = sampleMeanRay();
    const float_vec n = meanRay.d;
    const float_type theta = std::asin(2 * PCG32::rand() - 1);
    const float_type cos0 = std::cos(theta);
    const float_type sin0 = std::sin(theta);
    const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
    return std::pair(Ray(meanRay.o, rotDir), theta * RAD_TO_DEG);
}

std::pair<std::vector<Ray>, std::vector<Ray>> Circle::generateExtremeDiffuseRays() const {
    std::vector<Ray> p1ExtrRays, p2ExtrRays;
    for (int degrees = 0; degrees < 360; degrees += 2) {
        const float_type theta = degrees * DEG_TO_RAD;
        const float_vec p = c + r * float_vec(std::cos(theta), std::sin(theta));
        const float_vec n = (p - c) / r;
        const float_vec p1RotDir(-n.y, n.x);
        const float_vec p2RotDir(n.y, -n.x);
        p1ExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), p1RotDir);
        p2ExtrRays.emplace_back(p + (DOINKING ? DOINK * n : float_vec(0, 0)), p2RotDir);
    }
    return std::pair(p1ExtrRays, p2ExtrRays);
}

std::pair<std::vector<Ray>, std::vector<Ray>> Circle::generateExtremeInfiniteRays(const int &numRays) const {
    return std::pair(std::vector<Ray>(), std::vector<Ray>());
}

std::vector<Ray> Circle::generateFinalPlotRays() const {
    std::vector<Ray> rays;
    for (int degrees = -170; degrees <= 170; degrees += 20) {
        const float_type theta = degrees * DEG_TO_RAD;
        const float_type cos0 = std::cos(theta);
        const float_type sin0 = std::sin(theta);
        const float_vec n = float_vec(-1, 0);
        const float_vec rotDir(n.x * cos0 - n.y * sin0, n.x * sin0 + n.y * cos0);
        rays.emplace_back(getCentre() + ((DOINKING ? DOINK : 0) + r) * rotDir, rotDir);
    }
    return rays;
}
