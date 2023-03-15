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
    inline void init(size_t node_number, size_t bond_number, size_t part_number) {
      node_positions_.init(0, node_number);
      node_velocities_.init(0, node_number);
      node_adjacents_bonds1_.init(0, node_number);
      node_adjacents_bonds2_.init(0, node_number);
      node_adjacents_curvature_.init(0, node_number);
      part_node_emphasis_.init(0, part_number);
      part_node_rigid1_.init(0, part_number);
      part_node_rigid2_.init(0, part_number);
      bond_relations1_.init(0, bond_number);
      bond_relations2_.init(0, bond_number);
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
    inline CoreMath::Array<CoreMath::Vector> d_get_positions(int center, 
        const CoreMath::Array<int>& others) {
      CoreMath::Array<CoreMath::Vector> result(others.size());
      for (int i=0; i<others.size(); i++)
        result[i] = node_positions_(others[i])-node_positions_(center);
      return result;
    }
    inline CoreMath::Pair<CoreMath::Vector> d_get_positions(int center, 
        const CoreMath::Pair<int>& others) {
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
    inline void sync() {
      node_positions_.sync<Device>();
      node_velocities_.sync<Device>();
      node_adjacents_bonds1_.sync<Device>();
      node_adjacents_bonds2_.sync<Device>();
      node_adjacents_curvature_.sync<Device>();
      part_node_emphasis_.sync<Device>();
      part_node_rigid1_.sync<Device>();
      part_node_rigid2_.sync<Device>();
      bond_relations1_.sync<Device>();
      bond_relations2_.sync<Device>();
    }
    
    /// @brief alias of default memory and host mirrorspace
    using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
    using HostMirrorSpace = Kokkos::Impl::HostMirror<MemorySpace>::Space;

    void dump(std::string file_name, Metadata::DumpType dump_type);
    void store(std::string file_name);
    void load(std::string file_name);
    void update();

    /// @brief Random number pool
    CoreMath::Pool rand_pool_;

  // Data which will be store and load
  public:
    /// @brief These are used in calculate
    CoreMath::View<CoreMath::Vector> node_positions_;  
    CoreMath::View<CoreMath::Vector> node_velocities_;

    /// @brief Adjacents of nodes
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds1_;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_bonds2_;
    CoreMath::View<CoreMath::Array<int>> node_adjacents_curvature_;

    /// @brief Nodes of dislocations, have different types when output.
    CoreMath::View<int> part_node_emphasis_;

    /// @brief Nodes regards like rigid body, for edge processing.
    CoreMath::View<int> part_node_rigid1_;
    CoreMath::View<int> part_node_rigid2_;

    /// @brief These are used just for output, bonds and bond types
    CoreMath::View<CoreMath::Pair<int>> bond_relations1_;
    CoreMath::View<CoreMath::Pair<int>> bond_relations2_;
  
  private:
    int __time_step = 0;
}; // class ModelSystem

#endif // QUADRATUBE_MODEL_SYSTEM_H_
