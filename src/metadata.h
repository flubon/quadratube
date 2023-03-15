/**
 * @file metadata.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief metadata for dump and energy function
 * @version 0.0.1
 * @date 2023-03-14
 * 
 * @copyright Copyright (c) 2023
 */
#ifndef QUADRATUBE_METADATA_H_
#define QUADRATUBE_METADATA_H_

#include <cstdint>

#include "core/math.h"
#include "core/energy.h"

namespace Metadata {

/**
 * @brief meta data for energy funtions
 * @details Decide which function to use for update
 */
class EnergyMetaData {
  public:
    /// @brief alias for runtime change, don't use these names in class System,
    ///     in-class initializing can only use in c++11 (and after)
    double& bond1_rest_length_ = __data[0];
    double& bond1_spring_constant_ = __data[1];
    double& bond2_rest_length_ = __data[3];
    double& bond2_spring_constant_ = __data[4];
    double& bond2_truncate_length_ = __data[5];
    double& curvature_bending_rigidity_ = __data[6];

    /// @brief step length, mass of each particle, damping coefficient, temperature
    double step_length_ = 1e-2;
    double mass_ = 1;
    double damp_coeff_ = 1;
    double temperature_ = 0;

    /// @brief alias of energy function and gradient function for update and dump
    ///    function, use these in class System instead of direct CoreEnergy function
    inline double bond1_energy(const CoreMath::Vector& other) {
      return CoreEnergy::harmonic_energy(__data, other);
    }
    inline CoreMath::Vector bond1_gradient(const CoreMath::Vector& other) {
      return CoreEnergy::harmonic_gradient(__data, other);
    }
    inline double bond2_energy(const CoreMath::Vector& other) {
      return CoreEnergy::ljts_energy(__data+3, other);
    }
    inline CoreMath::Vector bond2_gradient(const CoreMath::Vector& other) {
      return CoreEnergy::ljts_gradient(__data+3, other);
    }
    inline double curvature_energy(const CoreMath::Array<CoreMath::Vector>& others) {
      return CoreEnergy::curvature_energy(__data+6, others);
    }
    inline CoreMath::Array<CoreMath::Vector> curvature_gradient(
        const CoreMath::Array<CoreMath::Vector>& others) {
      return CoreEnergy::curvature_gradient(__data+6, others);
    }
  
  private:
    /// @brief default parameter settings, we just need to change bond2_spring_constant
    ///     and curvature_bending_rigidity
    double __data[8] = {1, Kokkos::sqrt(3)/2, 0, Kokkos::sqrt(2), 0.1, 2.5*Kokkos::sqrt(2), 0.1};
};

typedef uint64_t DumpType;

/// @brief output all bonds into datafile (predifined)
const DumpType kPrintBond         = 1 << 0;
/// @brief don't output bonds of type 2 (predifined)
const DumpType kExcludeBondType2  = 1 << 1;
/// @brief give dislocation atoms a different type (predifined)
const DumpType kPrintDislocations = 1 << 2;
/// @brief output velocities of each atom (predifined)
const DumpType kPrintVelocities   = 1 << 3;

/// @brief custom marcos (defined by yourself)
const DumpType kPrintPotentialEnergy   = 1 << 4;
const DumpType kPrintGaussianCurvature = 1 << 5;
const DumpType kPrintMeanCurvature     = 1 << 6;

/// @brief alias of print all informations, avoid redundant macros
const DumpType kPrintAll = kPrintBond | kPrintDislocations |
    kPrintVelocities | kPrintPotentialEnergy |
    kPrintGaussianCurvature | kPrintMeanCurvature;

/// @brief simplified writing
typedef const CoreMath::Array<CoreMath::Vector>& ConstAdjacentNodes;

/**
 * @brief metadata for function dump compute (custom compute)
 */
struct {
  const char* name;
  DumpType dump_type;
  double (*func_num)(EnergyMetaData, ConstAdjacentNodes,
      ConstAdjacentNodes, ConstAdjacentNodes);
} const kDumpMetaData[] = {
  {"c_epot", kPrintPotentialEnergy, [](EnergyMetaData para, ConstAdjacentNodes others1,
      ConstAdjacentNodes others2, ConstAdjacentNodes others3){
    double count = 0;
    for (int i=0; i<others1.size(); i++)
      count += para.bond1_energy(others1[i]) / 2;
    for (int i=0; i<others2.size(); i++)
      count += para.bond2_energy(others2[i]) / 2;
    return para.curvature_energy(others3) + count;
  }},
  {"g_curv", kPrintGaussianCurvature, [](EnergyMetaData, ConstAdjacentNodes,
      ConstAdjacentNodes, ConstAdjacentNodes others){
    return CoreEnergy::gaussian_curvature(others);
  }},
  {"m_curv", kPrintMeanCurvature, [](EnergyMetaData, ConstAdjacentNodes,
      ConstAdjacentNodes, ConstAdjacentNodes others){
    return CoreEnergy::mean_curvature(others);
  }}
};

/// check whether name is in i
#define DUMP_CHECK(name, i) ((name & i) != 0)

/// 0: gradient descent, 1: langevin dynamics, 2: overdamped langevin dynamics
#define DYNAMICS 2

} // namespace Metadata

#endif // QUADRATUBE_METADATA_H_