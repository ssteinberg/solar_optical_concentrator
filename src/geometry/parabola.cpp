#include "geometry/parabola.h"

float_vec Parabola::getCentre() const { return float_vec(c, 0); }

float_type Parabola::getLength() const { return float_type(0); }

bool Parabola::intersect(const Ray &ray, HitInfo &hitInfo) const {
    if (IGNORE_SOURCE && type == Type::SOURCE) return false;
    else {
        if (ray.d.y == 0) {
            const float_type s = (a * ray.o.y * ray.o.y - ray.o.x + c) / ray.d.x;
            if (s < 0) return false;
            hitInfo.l = s;
            hitInfo.p = ray.o + hitInfo.l * ray.d;
            hitInfo.n = normalize(float_vec(-1, 2 * a * hitInfo.p.y));
            hitInfo.t = type;
            return true;
        }
        const float_type A = a * ray.d.y * ray.d.y;
        const float_type B = 2 * a * ray.o.y * ray.d.y - ray.d.x;
        const float_type C = a * ray.o.y * ray.o.y + c - ray.o.x;
        const float_type det = B * B - 4 * A * C;
        if (det < 0) return false;
        const float_type s1 = (-B - std::sqrt(det)) / (2 * A);
        const float_type s2 = (-B + std::sqrt(det)) / (2 * A);
        float_type s = 0;
        if (s1 >= 0) s = s1;
        else if (s2 >= 0) s = s2;
        else return false;
        hitInfo.l = s;
        hitInfo.p = ray.o + hitInfo.l * ray.d;
        hitInfo.n = normalize(float_vec(-1, 2 * a * hitInfo.p.y));
        hitInfo.t = type;
        return true;
    }
}

Ray Parabola::sampleMeanRay() const {
    return Ray(float_vec(0, 0), float_vec(0, 0));
}

std::pair<Ray, float_type> Parabola::sampleDiffuseRay() const {
    return std::pair(Ray(float_vec(0, 0), float_vec(0, 0)), float_type(0));
}

std::pair<std::vector<Ray>, std::vector<Ray>> Parabola::generateExtremeDiffuseRays() const {
    return std::pair(std::vector<Ray>(0), std::vector<Ray>(0));
}

std::pair<std::vector<Ray>, std::vector<Ray>> Parabola::generateExtremeInfiniteRays(const int &numRays) const {
    return std::pair(std::vector<Ray>(0), std::vector<Ray>(0));
}

std::vector<Ray> Parabola::generateFinalPlotRays() const {
    return std::vector<Ray>();
}
