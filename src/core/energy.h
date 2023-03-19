/**
 * @file energy.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-08
 * 
 * @copyright Copyright (c) 2023
 */
#ifndef QUADRATUBE_CORE_ENERGY_H_
#define QUADRATUBE_CORE_ENERGY_H_

#include <Kokkos_Core.hpp>

#include "core/math.h"

namespace CoreEnergy {

/**
 * @brief harmonic bond energy
 * 
 * @param para para[0]: rest length, para[1]: elastic coefficient k
 * @param other 
 * @return double 
 */
KOKKOS_INLINE_FUNCTION
double harmonic_energy(const double* para, const CoreMath::Vector& other) {
  return para[1]/2 * (CoreMath::mod(other)-para[0]) * (CoreMath::mod(other)-para[0]);
}
KOKKOS_INLINE_FUNCTION
CoreMath::Vector harmonic_gradient(const double* para, const CoreMath::Vector& other) {
  return para[1] * (1-para[0]/CoreMath::mod(other)) * other;
}

/**
 * @brief Lennard-Jones truncated & shifted potential
 * 
 * @param para para[0]: rest length, para[1]: coefficient, para[2]: cutoff length
 * @param other 
 * @return double 
 */
KOKKOS_INLINE_FUNCTION
double ljts_energy(const double* para, const CoreMath::Vector& other) {
  if (CoreMath::mod(other) > para[2])
    return 0;
  // (r_min / r_end)^6
  double rmin_rend6 = Kokkos::pow(para[0]/para[2], 6);
  double cutoff_energy = (rmin_rend6 - 2) * rmin_rend6;
  // (r_min / r)^6
  double rmin_r6 = Kokkos::pow(para[0]/CoreMath::mod(other), 6);
  return para[1] * (rmin_r6 - 2) * rmin_r6 - cutoff_energy;
}
KOKKOS_INLINE_FUNCTION
CoreMath::Vector ljts_gradient(const double* para, const CoreMath::Vector& other) {
  if (CoreMath::mod(other) > para[2])
    return CoreMath::Vector();
  // r_min^6 / r^7
  double rmin6_r7 = Kokkos::pow(para[0], 6) / Kokkos::pow(CoreMath::mod(other), 7);
  return 12 * para[1] * rmin6_r7*(1/CoreMath::mod(other) - rmin6_r7) * other;
}

/**
 * @brief mean curvature and gaussian curvature
 * 
 * @param others 
 * @return double 
 */
KOKKOS_FUNCTION
double mean_curvature(const CoreMath::Array<CoreMath::Vector>& others);
KOKKOS_FUNCTION
double gaussian_curvature(const CoreMath::Array<CoreMath::Vector>& others);

/**
 * @brief total energy arised from curvature
 * 
 * @param para para[0]: coefficient
 * @param others 
 * @return double 
 */
KOKKOS_FUNCTION
double curvature_energy(const double* para, const CoreMath::Array<CoreMath::Vector>& others);
KOKKOS_FUNCTION
CoreMath::Array<CoreMath::Vector> curvature_gradient(const double* para, 
    const CoreMath::Array<CoreMath::Vector>& others);

} // namespace CoreEnergy

#endif // QUADRATUBE_CORE_ENERGY_H_