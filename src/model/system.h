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
    /// @brief initialize members
    /// @param node_number 
    /// @param bond_number 
    /// @param part_number 
    inline void init_all(size_t node_number, size_t bond_number) {
      for (int i=0; i<2; i++)
        __node_attributes[i].init(0, node_number);
      for (int i=0; i<3; i++)
        __node_relations[i].init(0, node_number);
      for (int i=0; i<5; i++)
        __node_judgements[i].init(0, node_number);
      for (int i=0; i<2; i++)
        __bond_relations[i].init(0, bond_number);
    }

    /// @brief make nodes number consistent
    /// @param node_number 
    inline void resize_nodes(size_t node_number) {
      for (int i=0; i<2; i++)
        __node_attributes[i].resize(node_number);
      for (int i=0; i<3; i++)
        __node_relations[i].resize(node_number);
      for (int i=0; i<5; i++)
        __node_judgements[i].resize(node_number);
    }

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
    
    /**
     * @brief sync of all Kokkos::DualView
     * 
     * @tparam Device 
     */
    template <class Device>
    inline void sync_all() {
      for (int i=0; i<2; i++)
        __node_attributes[i].sync<Device>();
      for (int i=0; i<3; i++)
        __node_relations[i].sync<Device>();
      for (int i=0; i<5; i++)
        __node_judgements[i].sync<Device>();
      for (int i=0; i<2; i++)
        __bond_relations[i].sync<Device>();
    }

    template <class Device>
    inline void modify_all() {
      for (int i=0; i<2; i++)
        __node_attributes[i].modify<Device>();
      for (int i=0; i<3; i++)
        __node_relations[i].modify<Device>();
      for (int i=0; i<5; i++)
        __node_judgements[i].modify<Device>();
      for (int i=0; i<2; i++)
        __bond_relations[i].modify<Device>();
    }
    
    /// @brief alias of default memory and host mirrorspace
    using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
    using HostMirrorSpace = Kokkos::Impl::HostMirror<MemorySpace>::Space;

    void dump(std::string file_name, Metadata::DumpType dump_type);
    void store(std::string file_name);
    void load(std::string file_name);
    void update(bool just_velocity = false);

    /// @brief Random number pool
    CoreMath::Pool rand_pool_;

  // Data which will be store and load
  public:
    /// @brief These are used in calculate
    CoreMath::View<CoreMath::Vector>& node_positions_ = __node_attributes[0];  
    CoreMath::View<CoreMath::Vector>& node_velocities_ = __node_attributes[1];

    /// @brief Adjacents of nodes
    CoreMath::View<CoreMath::Array<int>>& node_adjacents_bonds1_ = __node_relations[0];
    CoreMath::View<CoreMath::Array<int>>& node_adjacents_bonds2_ = __node_relations[1];
    CoreMath::View<CoreMath::Array<int>>& node_adjacents_curvature_ = __node_relations[2];

    /// @brief Nodes of dislocations, have different types when output.
    CoreMath::View<bool>& node_if_emphasis_ = __node_judgements[0];
    int& node_if_emphasis_count_ = __node_judgements_count[0];

    /// @brief Nodes regards like rigid body, for edge processing.
    CoreMath::View<bool>& node_if_rigid1_ = __node_judgements[1];
    int& node_if_rigid1_count_ = __node_judgements_count[1];
    CoreMath::View<bool>& node_if_next_to_rigid1_ = __node_judgements[2];
    int& node_if_next_to_rigid1_count_ = __node_judgements_count[2];
    CoreMath::View<bool>& node_if_rigid2_ = __node_judgements[3];
    int& node_if_rigid2_count_ = __node_judgements_count[3];
    CoreMath::View<bool>& node_if_next_to_rigid2_ = __node_judgements[4];
    int& node_if_next_to_rigid2_count_ = __node_judgements_count[4];

    /// @brief These are used just for output, bonds and bond types
    CoreMath::View<CoreMath::Pair<int>>& bond_relations1_ = __bond_relations[0];
    CoreMath::View<CoreMath::Pair<int>>& bond_relations2_ = __bond_relations[1];
  
  private:
    CoreMath::View<CoreMath::Vector> __node_attributes[2];
    CoreMath::View<CoreMath::Array<int>> __node_relations[3];
    CoreMath::View<bool> __node_judgements[5];
    CoreMath::View<CoreMath::Pair<int>> __bond_relations[2];
    /// @brief count how many trues in __node_judgements
    int __node_judgements_count[5] = {0};

    int __time_step = 0;
}; // class ModelSystem

#endif // QUADRATUBE_MODEL_SYSTEM_H_
