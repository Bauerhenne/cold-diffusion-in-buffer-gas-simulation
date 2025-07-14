#include "flowfieldentries.h"

bool FlowFieldEntries::operator==(FlowFieldEntries const &rhs)
{
    if(v[0]==rhs.v[0] && v[1]==rhs.v[1] && v[2]==rhs.v[2] && rho==rhs.rho && T==rhs.T)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FlowFieldEntries::operator!=(FlowFieldEntries const &rhs)
{
    if(*this == rhs)
    {
        return false;
    }
    else
    {
        return true;
    }
}

FlowFieldEntries& FlowFieldEntries::operator=(FlowFieldEntries const &rhs)
{
    v[0] = rhs.v[0];
    v[1] = rhs.v[1];
    v[2] = rhs.v[2];
    rho  = rhs.rho;
    T    = rhs.T;
    return *this;
}

FlowFieldEntries& FlowFieldEntries::operator+=(FlowFieldEntries const &rhs)
{
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    rho  += rhs.rho;
    T    += rhs.T;
    return *this;
}

FlowFieldEntries& FlowFieldEntries::operator-=(FlowFieldEntries const &rhs)
{
    v[0] -= rhs.v[0];
    v[1] -= rhs.v[1];
    v[2] -= rhs.v[2];
    rho  -= rhs.rho;
    T    -= rhs.T;
    return *this;
}

FlowFieldEntries& FlowFieldEntries::operator*=(double const &rhs)
{
    v[0] *= rhs;
    v[1] *= rhs;
    v[2] *= rhs;
    rho  *= rhs;
    T    *= rhs;
    return *this;
}

std::ostream& operator<<(std::ostream &stream, FlowFieldEntries const &rhs)
{
    return stream << rhs.v[0] << " " << rhs.v[1] << " " << rhs.v[2] << " " << rhs.rho << " " << rhs.T << "\n";
}

std::istream& operator>>(std::istream &stream, FlowFieldEntries &rhs)
{
    return stream >> rhs.v[0] >> rhs.v[1] >> rhs.v[2] >> rhs.rho >> rhs.T;
}

void FlowFieldEntries::set_zero()
{
    v[0] = 0;
    v[2] = 0;
    v[1] = 0;
    rho = 0;
    T = 0;
}
