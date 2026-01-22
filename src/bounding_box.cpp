#include "bounding_box.h"

BoundingBox::BoundingBox(): min(std::numeric_limits<float_type>::min(), std::numeric_limits<float_type>::min()),
                            max(std::numeric_limits<float_type>::max(), std::numeric_limits<float_type>::max()) {}

void BoundingBox::fit(const LineSegment &l) {
    if (l.p1.x < min.x) min.x = l.p1.x;
    if (l.p1.y < min.y) min.y = l.p1.y;
    if (l.p2.x < min.x) min.x = l.p2.x;
    if (l.p2.y < min.y) min.y = l.p2.y;
    if (l.p1.x > max.x) max.x = l.p1.x;
    if (l.p1.y > max.y) max.y = l.p1.y;
    if (l.p2.x > max.x) max.x = l.p2.x;
    if (l.p2.y > max.y) max.y = l.p2.y;
}

bool BoundingBox::intersect(const Ray &ray, HitInfo &minHit) const {
    float_type tx1 = (min.x - ray.o.x) / ray.d.x;
    float_type ty1 = (min.y - ray.o.y) / ray.d.y;
    float_type tx2 = (max.x - ray.o.x) / ray.d.x;
    float_type ty2 = (max.y - ray.o.y) / ray.d.y;
    if (tx1 > tx2) {
        const float_type temp = tx1;
        tx1 = tx2;
        tx2 = temp;
    }
    if (ty1 > ty2) {
        const float_type temp = ty1;
        ty1 = ty2;
        ty2 = temp;
    }
    float_type t1 = tx1;
    if (t1 < ty1) t1 = ty1;
    float_type t2 = tx2;
    if (t2 > ty2) t2 = ty2;
    if (t1 > t2) return false;
    if ((t1 < 0) && (t2 < 0)) return false;
    minHit.l = t1;
    return true;
}
