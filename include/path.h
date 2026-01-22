#pragma once

#include "types.h"

struct Path {
    std::vector<float_vec> vertices;

    void addVertex(const float_vec &v);

    void writePath(std::ofstream &file) const;

    void writeFinalSegment(std::ofstream &file) const;
};
