#include <fstream>

#include "path.h"
#include "types.h"

void Path::addVertex(const float_vec &v) {
    vertices.push_back(v);
}

void Path::writePath(std::ofstream &file) const {
    const int pathLength = vertices.size() - 1;
    if (pathLength > 0) {
        for (int i = 0; i < pathLength; ++i) {
            const float_vec v1 = vertices[i];
            const float_vec v2 = vertices[i + 1];
            file << v1.x << "," << v1.y << "," << v2.x << "," << v2.y << "\n";
        }
    }
}

void Path::writeFinalSegment(std::ofstream &file) const {
    const int pathLength = vertices.size() - 1;
    if (pathLength > 0) {
        const float_vec v1 = vertices[pathLength - 1];
        const float_vec v2 = vertices[pathLength];
        file << v1.x << "," << v1.y << "," << v2.x << "," << v2.y << "\n";
    }
}
