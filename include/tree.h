#pragma once

#include "node.h"

struct Tree {
    Node* root = nullptr;
    Tree() = default;
    Tree(const std::vector<LineSegment>& v);

    bool intersect(const Ray& ray, HitInfo& hitInfo) const;
};
