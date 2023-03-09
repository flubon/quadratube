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

#include "core/math.h"

namespace CoreEnergy {

/**
 * @brief local harmonic bond energy, 1/2 of all bonds' energy around the node
 * 
 * @param para para[0]: rest lenth, para[1]: elastic coefficient k
 * @param center 
 * @param others 
 * @return double 
 */
double harmonic_bond_energy(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief local force arised from harmonic bond energy
 * 
 * @param para para[0]: rest lenth, para[1]: elastic coefficient k
 * @param center 
 * @param others 
 * @return CoreMath::Vector 
 */
CoreMath::Vector harmonic_bond_force(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief local Lennard-Jones truncated & shifted potential
 * 
 * @param para para[0]: rest lenth, para[1]: coefficient, para[2]: cutoff lenth
 * @param center 
 * @param others 
 * @return double 
 */
double ljts_bond_energy(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief local force arised from LJTS potential
 * 
 * @param para para[0]: rest lenth, para[1]: coefficient, para[2]: cutoff lenth
 * @param center 
 * @param others 
 * @return CoreMath::Vector 
 */
CoreMath::Vector ljts_bond_force(CoreMath::Array<double> para, CoreMath::Vector center, 
    CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief 
 * 
 * @param center 
 * @param others 
 * @return double 
 */
double mean_node_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief 
 * 
 * @param center 
 * @param others 
 * @return double 
 */
double gaussian_node_curvature(CoreMath::Vector center, CoreMath::Array<CoreMath::Vector> others);

/**
 * @brief 
 * 
 * @param para 
 * @param center 
 * @param others 
 * @return double 
 */
double curvature_node_energy(CoreMath::Array<double> para, CoreMath::Vector center,
    CoreMath::Array<CoreMath::Vector> others);

} // namespace CoreEnergy

#endif // QUADRATUBE_CORE_ENERGY_H_