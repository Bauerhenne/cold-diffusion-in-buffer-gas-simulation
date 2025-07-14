#pragma once
#include <iostream>

class Vec3d {
public:
    double r[3]{ 0,0,0 };

    bool operator==(Vec3d const& rhs);
    bool operator!=(Vec3d const& rhs);

    friend std::ostream& operator<<(std::ostream& stream, Vec3d const& rhs);
    friend std::istream& operator>>(std::istream& stream, Vec3d& rhs);
};