#include "tree.h"

Tree::Tree(const std::vector<LineSegment> &v) { root = new Node(v, 0, static_cast<int>(v.size()) - 1); }

bool Tree::intersect(const Ray &ray, HitInfo &hitInfo) const { if (root) return root->intersect(ray, hitInfo); return false; }
