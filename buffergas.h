#pragma once
#include "constants.h"
#include "flowfieldentries.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include <math.h>

class BufferGas
{
public:
    double mass;        // unit kg
    double temperature; // unit K
    double density;     // unit particle/m^3
    double vbThermSquare;// average thermal velocity square from Maxwell Boltzmann distribution
    double vDrift[3];

    BufferGas();

    void load_Data(const std::string filename);
    bool is_in_volume(double const coord[3]);
    bool hits_exit(double const coord[3]);
    void update_coord(const double* coord);
    void get_thermal_velocity(double* u);
    double get_radius();

private:
    std::vector<FlowFieldEntries> fField;
    int numData = 0;
    int numDataEfficient = 0;
    int nx = 0, ny = 0, nz = 0;
    int ixold = -1, iyold = -1, izold = -1;
    double delta = 0.0005;
    double diameter = 0;
    double exitDiameter = 0;
    double halfDiameter = 0;
    double halfDiameterSquare;
    double height;
    std::default_random_engine generator;
    std::normal_distribution<double> normalDistribution;
};