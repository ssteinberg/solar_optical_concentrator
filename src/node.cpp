#include "node.h"

Node::Node(const std::vector<LineSegment> &v, const int &li, const int &ri): size(ri - li + 1) {
    if (size > 1) {
        for (int i = li; i <= ri; ++i) bb.fit(v[i]);
        left = new Node(v, li, li + size / 2 - 1);
        right = new Node(v, li + size / 2, ri);
    }
    else ls = &v.at(li);
}

bool Node::intersect(const Ray &ray, HitInfo &hitInfo) const {
    if (size == 1) return ls->intersect(ray, hitInfo);
    if (bb.intersect(ray, hitInfo)) {
        if (left->intersect(ray, hitInfo)) return true;
        if (right->intersect(ray, hitInfo)) return true;
    }
    return false;
}
