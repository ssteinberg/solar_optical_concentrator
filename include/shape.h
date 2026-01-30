#pragma once

#include <iostream>

enum struct Shape {
    FLAT,
    CYLINDRICAL,
    ELLIPTICAL,
    PARABOLIC,
    CONSTRUCTED
};

inline std::ostream& operator<<(std::ostream& os, const Shape shape)
{
    return os << static_cast<int>(shape);
}
