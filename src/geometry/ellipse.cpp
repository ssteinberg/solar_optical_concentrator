#include "geometry/ellipse.h"

float_vec Ellipse::getCentre() const { return float_vec(-L / 2, 0); }

float_type Ellipse::getLength() const { return float_type(0); }

bool Ellipse::intersect(const Ray &ray, HitInfo &hitInfo) const {
    if (IGNORE_SOURCE && type == Type::SOURCE) return false;
    else {
        const float_vec centre = float_vec(-L / 2, 0);
        const float_vec oc = ray.o - centre;
        float_type inv_a2 = 1 / (saa * saa);
        float_type inv_b2 = 1 / (sab * sab);
        const float_type a = ray.d.x * ray.d.x * inv_a2 + ray.d.y * ray.d.y * inv_b2;
        const float_type b = 2 * (oc.x * ray.d.x * inv_a2 + oc.y * ray.d.y * inv_b2);
        const float_type c = oc.x * oc.x * inv_a2 + oc.y * oc.y * inv_b2 - 1;
        const float_type d = b * b - 4 * a * c;
        if (d < 0) return false;
        const float_type s1 = (-b - std::sqrt(d)) / (2 * a);
        const float_type s2 = (-b + std::sqrt(d)) / (2 * a);
        float_type s = 0;
        if (s1 >= 0) s = s1;
        else if (s2 >= 0) s = s2;
        else return false;
        hitInfo.l = s;
        hitInfo.p = ray.o + hitInfo.l * ray.d;
        hitInfo.n = -normalize(float_vec(hitInfo.p.x - (-L / 2), (saa * saa / (sab * sab)) * hitInfo.p.y));
        hitInfo.t = type;
        return true;
    }
}

Ray Ellipse::sampleMeanRay() const {
    return Ray(float_vec(0, 0), float_vec(0, 0));
}

std::pair<Ray, float_type> Ellipse::sampleDiffuseRay() const {
    return std::pair(Ray(float_vec(0, 0), float_vec(0, 0)), float_type(0));
}

std::pair<std::vector<Ray>, std::vector<Ray>> Ellipse::generateExtremeDiffuseRays() const {
    return std::pair(std::vector<Ray>(0), std::vector<Ray>(0));
}

std::pair<std::vector<Ray>, std::vector<Ray>> Ellipse::generateExtremeInfiniteRays(const int &numRays) const {
    return std::pair(std::vector<Ray>(0), std::vector<Ray>(0));
}

std::vector<Ray> Ellipse::generateFinalPlotRays() const {
    return std::vector<Ray>();
}
