#include "vec3d.h"

bool Vec3d::operator==(Vec3d const& rhs)
{
    if (r[0] == rhs.r[0] && r[1] == rhs.r[1] && r[2] == rhs.r[2])
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Vec3d::operator!=(Vec3d const& rhs)
{
    if (*this == rhs)
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::ostream& operator<<(std::ostream& stream, Vec3d const& rhs)
{
    return stream << rhs.r[0] << " " << rhs.r[1] << " " << rhs.r[2];
}

std::istream& operator>>(std::istream& stream, Vec3d& rhs)
{
    return stream >> rhs.r[0] >> rhs.r[1] >> rhs.r[2];
}
