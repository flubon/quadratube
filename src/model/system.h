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

#include <string>

#include "metadata.h"
#include "core/math.h"

/**
 * @class System
 * @brief Molecular system
 * @details The smallest model of molecular.
 * 
 */
class ModelSystem : public Metadata::EnergyMetaData {
  public:
    /// @brief transform index of others to relative positions
    /// @details "h" is for host, "d" is for device
    /// @param center 
    /// @param others 
    /// @return
    inline CoreMath::Array<CoreMath::Vector> h_get_positions(int center, 
        const CoreMath::Array<int>& others) {
      CoreMath::Array<CoreMath::Vector> result(others.size());
      for (int i=0; i<others.size(); i++)
        result[i] = node_positions_[others[i]]-node_positions_[center];
      return result;
    }
    inline CoreMath::Pair<CoreMath::Vector> h_get_positions(int center, 
        const CoreMath::Pair<int>& others) {
      return CoreMath::Pair<CoreMath::Vector>(
        node_positions_[others[0]]-node_positions_[center], 
        node_positions_[others[1]]-node_positions_[center]);
    }
    KOKKOS_INLINE_FUNCTION
    CoreMath::Array<CoreMath::Vector> d_get_positions(int center, 
        const CoreMath::Array<int>& others) const {
      CoreMath::Array<CoreMath::Vector> result(others.size());
      for (int i=0; i<others.size(); i++)
        result[i] = node_positions_(others[i])-node_positions_(center);
      return result;
    }
    KOKKOS_INLINE_FUNCTION
    CoreMath::Pair<CoreMath::Vector> d_get_positions(int center, 
        const CoreMath::Pair<int>& others) const {
      return CoreMath::Pair<CoreMath::Vector>(
        node_positions_(others[0])-node_positions_(center), 
        node_positions_(others[1])-node_positions_(center));
    }
    
    /// @brief alias of default memory and host mirrorspace
    using MemorySpace = CoreMath::View<int>::MemorySpace;
    using HostMirrorSpace = CoreMath::View<int>::HostMirrorSpace;

    void dump(std::string file_name, Metadata::DumpType dump_type);
    void store(std::string file_name);
    void load(std::string file_name);
    void update(bool just_velocity = false);

    /// @brief Random number pool
    CoreMath::Pool rand_pool_;

  // Data which will be store and load
  public:
    /// @brief These are used in calculate
    CoreMath::View<CoreMath::Vector> node_positions_;  
    CoreMath::View<CoreMath::Vector> node_velocities_;

    // These won't be changed during update, and won't be checked by update or store.

    /// @brief Adjacents of nodes
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds1_;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds2_;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_curvature_;

    /// @brief Nodes of dislocations, have different types when output.
    CoreMath::View<bool> node_if_emphasis_;
    int node_if_emphasis_count_ = 0;

    /// @brief Nodes regards like rigid body, for edge processing.
    CoreMath::View<bool> node_if_rigid1_;
    int node_if_rigid1_count_ = 0;
    CoreMath::View<bool> node_if_next_to_rigid1_;
    int node_if_next_to_rigid1_count_ = 0;
    CoreMath::View<bool> node_if_rigid2_;
    int node_if_rigid2_count_ = 0;
    CoreMath::View<bool> node_if_next_to_rigid2_;
    int node_if_next_to_rigid2_count_ = 0;

    /// @brief These are used just for output, bonds and bond types
    CoreMath::View<CoreMath::Pair<int>> bond_relations1_;
    CoreMath::View<CoreMath::Pair<int>> bond_relations2_;
  
  private:
    int __time_step = 0;
}; // class ModelSystem

#endif // QUADRATUBE_MODEL_SYSTEM_H_
