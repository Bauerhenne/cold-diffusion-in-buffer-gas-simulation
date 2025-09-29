#pragma once
#include "constants.h"
#include <iostream>
#include <time.h>
#include <random>
#include <math.h>
#include <fstream>

class Molecule
{
 public:
  double coord[3], velocity[3];
  double mass; //unit kg

  Molecule();

  void init_coordinates_velocity(const double temperature, double radius);
  double get_norm_vrel(double vb[3]);
  double get_temperature(double vbDrift[3]);
  double get_uniform_distributed_variable();
  void set_velocity_from_hard_sphere_collision_with_buffer_gas(double vb[3], double bufferGasMass);

 private:
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniformDistribution;
  std::normal_distribution<double> normalDistribution;
  double inverse_cos_square_dist(const double u);
  void check_cos_square_dist();
};
