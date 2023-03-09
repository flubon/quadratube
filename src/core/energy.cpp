/**
 * @file energy.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-08
 * 
 * @copyright Copyright (c) 2023
 */
#include "core/energy.h"

#include "core/math.h"

namespace CoreEnergy {

namespace {

double __mean_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others,
    double size) {

  CoreMath::Vector H_vec;

  for (int i = 0; i < others.size(); i++) {
    int left = (i==0) ? others.size()-1 : i-1;
    int right = (i==others.size()-1) ? 0 : i+1;

    // vectors for a node and its four related bonds
    CoreMath::Vector left_up = others[i] - others[left];
    CoreMath::Vector left_down = center - others[left];
    CoreMath::Vector right_up = others[i] - others[right];
    CoreMath::Vector right_down = center - others[right];

    // tan(acos(x)) = sqrt(1/x/x-1)
    double cos_angle1 = left_up*left_down/CoreMath::mod(left_up)/CoreMath::mod(left_down);
    double cos_angle2 = right_up*right_down/CoreMath::mod(right_up)/CoreMath::mod(right_down);
    H_vec += (1/Kokkos::sqrt(1/cos_angle1/cos_angle1-1) + 
      1/Kokkos::sqrt(1/cos_angle2/cos_angle2-1)) * (others[0] - center);
  }

  return CoreMath::mod(H_vec)/4/size;
}

double __gaussian_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others, 
    double size) {
  // Size around this node, total angle
  double angles = 2*PI;
  for (int i=0; i<others.size()-1; i++) {
    CoreMath::Vector v1 = others[i] - center;
    CoreMath::Vector v2 = others[i+1] - center;
    angles -= Kokkos::acos(v1*v2/CoreMath::mod(v1)/CoreMath::mod(v2));
  }
  CoreMath::Vector v1 = others[0] - center;
  CoreMath::Vector v2 = *(others.end()-1) - center;
  angles -= Kokkos::acos(v1*v2/CoreMath::mod(v1)/CoreMath::mod(v2));
  return angles/size;
}

double __size(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others) {
  double size = 0;
  for (int i=0; i<others.size()-1; i++)
    size += CoreMath::mod(CoreMath::cross(others[i]-center, others[i+1]-center)) / 2;
  size += CoreMath::mod(CoreMath::cross(others[0]-center, *(others.end()-1)-center)) / 2;
  return size/3;
}

} // namespace

double harmonic_bond_energy(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others) {
  double count = 0;
  for (auto i : others) {
    double r = CoreMath::mod(i-center);
    count += (r-para[0]) * (r-para[0]);
  }
  return para[1]/2 * count;
}

CoreMath::Vector harmonic_bond_force(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others) {
  CoreMath::Vector count;
  for (auto i : others) {
    CoreMath::Vector r = i - center;
    count += (1-para[0]/CoreMath::mod(r)) * r;
  }
  return para[1] * count;
}

double ljts_bond_energy(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others) {
  // (r_min / r_end)^6
  double rmin_rend6 = Kokkos::pow(para[0]/para[2], 6);
  double cutoff_energy = (rmin_rend6 - 2) * rmin_rend6;
  double count = 0;
  for (auto i : others) {
    double r = CoreMath::mod(i-center);
    double rmin_r6 = Kokkos::pow(para[0]/r, 6);
    count += (rmin_r6 - 2) * rmin_r6 - cutoff_energy;
  }
  return para[1] * count;
}

CoreMath::Vector ljts_bond_force(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others) {
  CoreMath::Vector count;
  for (auto i : others) {
    CoreMath::Vector r_vec = i-center;
    double r = CoreMath::mod(r_vec);
    // r_min^6 / r^7
    double rmin6_r7 = Kokkos::pow(para[0], 6) / Kokkos::pow(r, 7);
    count += rmin6_r7*(1/r - rmin6_r7) * r_vec;
  }
  return 12 * para[1] * count;
}

double mean_node_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others) {
  double size = __size(center, others);
  return __mean_curvature(center, others, size);
}

double gaussian_node_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others) {
  double size = __size(center, others);
  return __gaussian_curvature(center, others, size);
}
    
double curvature_node_energy(CoreMath::Array<double> para, CoreMath::Vector center,
    CoreMath::Array<CoreMath::Vector> others) {
  double size = __size(center, others);
  double mean = __mean_curvature(center, others, size);
  double gaussian = __gaussian_curvature(center, others, size);
  return para[0] * (2*mean*mean - gaussian);
}

} // namespace CoreEnergy