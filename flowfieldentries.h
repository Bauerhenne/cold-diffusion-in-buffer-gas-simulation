#pragma once
#include<iostream>

class FlowFieldEntries{
    public:

    double v[3]{0,0,0};
    double rho{0};
    double T{0};

    bool operator==(FlowFieldEntries const &rhs);
    bool operator!=(FlowFieldEntries const &rhs);
    FlowFieldEntries& operator=(FlowFieldEntries const &rhs);
    FlowFieldEntries& operator+=(FlowFieldEntries const &rhs);
    FlowFieldEntries& operator-=(FlowFieldEntries const &rhs);
    FlowFieldEntries& operator*=(double const &rhs);

    friend std::ostream& operator<<(std::ostream &stream, FlowFieldEntries const &rhs);
    friend std::istream& operator>>(std::istream &stream, FlowFieldEntries &rhs);

    void set_zero();
};