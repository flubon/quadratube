/**
 * @file system.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief Molecular system
 * @version 0.0.1
 * @date 2023-03-07
 * 
 * @copyright Copyright (c) 2023
 */
#ifndef QUADRATUBE_MODEL_SYSTEM_H_
#define QUADRATUBE_MODEL_SYSTEM_H_

#include "core/math.h"

namespace ModelSystem {

/**
 * @class System
 * @brief Molecular system
 * @details The smallest model of molecular.
 * 
 */
class System {
  public:
    void dump();
    void store();
    void load();

  public:
    /// @brief These are used in calculate
    CoreMath::View<CoreMath::Vector> node_positions;  
    CoreMath::View<CoreMath::Vector> node_velocities;

    /// @brief Adjacents of nodes
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds1;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds2;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_curvature;

    /// @brief Nodes of dislocations, have different types when output.
    CoreMath::View<int> part_node_emphasis;

    /// @brief Nodes regards like rigid body, for edge processing.
    CoreMath::View<int> part_node_rigid1;
    CoreMath::View<int> part_node_rigid2;

    /// @brief These are used just for output, bonds and bond types
    CoreMath::View<CoreMath::Pair<int>> bond_relations;
    CoreMath::View<int> bond_types;
};

} // namespace ModelSystem

#endif // QUADRATUBE_MODEL_SYSTEM_H_
