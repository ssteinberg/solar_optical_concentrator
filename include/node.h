#pragma once

#include "bounding_box.h"

struct Node {
    int size;
    BoundingBox bb;
    Node* left = nullptr;
    Node* right = nullptr;
    const LineSegment* ls = nullptr;

    Node(const std::vector<LineSegment>& v, const int& li, const int& ri);

    bool intersect(const Ray& ray, HitInfo& hitInfo) const;
};
