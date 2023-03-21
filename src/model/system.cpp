/**
 * @file system.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief Molecular system
 * @version 0.0.1
 * @date 2023-03-07
 * 
 * @copyright Copyright (c) 2023
 */
#include "model/system.h"

#include <stdio.h>

#include <Kokkos_Core.hpp>

#include "core/math.h"
#include "metadata.h"

namespace {

const char __data_file_header[] =
    "# Model for quadratube. AUTO generated, DO NOT EDIT\n\n"           // header
    "%li\tatoms\n%li\tbonds\n\n%i\tatom types\n%i\tbond types\n\n"        // atom & bond number
    "%.8f\t%.8f\txlo xhi\n%.8f\t%.8f\tylo yhi\n%.8f\t%.8f\tzlo zhi\n\n" // boundary
    "Masses\n\n1\t%.8f\n";

const char __dump_file_header[] = 
    "ITEM: TIMESTEP\n%i\nITEM: NUMBER OF ATOMS\n%li\nITEM: BOX BOUNDS ss ss ss\n"
    "%.8f\t%.8f\n%.8f\t%.8f\n%.8f\t%.8f\nITEM: ATOMS id type x y z";

} // namespace

void ModelSystem::dump(std::string file_name, Metadata::DumpType dump_type) {
  // others won't be changed by update
  node_positions_.sync<HostMirrorSpace>();  
  node_velocities_.sync<HostMirrorSpace>();

  // boundary: x, y, z: [0, 30*bond1_rest_length_]
  double boundary_max = 30*bond1_rest_length_;

  // total types. default type is 1, dislocations use type 2
  int atom_types = DUMP_CHECK(Metadata::kPrintDislocations, dump_type) ? 2 : 1;
  auto get_type = [=](bool in) {
    return (atom_types == 2) ? static_cast<int>(in) + 1 : 1;
  };

  // print bonds into data file, skip if file exists already.
  FILE *file = std::fopen((file_name + ".data").c_str(), "r");
  if (DUMP_CHECK(Metadata::kPrintBond, dump_type) || file == NULL) {
    file = std::fopen((file_name + ".data").c_str(), "w");

    // header, check whether bond_relations2_.size() == 0 is for model3
    size_t bonds = DUMP_CHECK(Metadata::kExcludeBondType2, dump_type) ?
        bond_relations1_.size() : bond_relations1_.size() + bond_relations2_.size();
    int bond_types = (bond_relations2_.size() == 0 ||
        DUMP_CHECK(Metadata::kExcludeBondType2, dump_type)) ? 1 : 2;
    std::fprintf(file, __data_file_header,
      node_positions_.size(), bonds, atom_types, bond_types,          // atom & bond number
      0., boundary_max, 0., boundary_max, 0., boundary_max, mass_  // boundary & masses
    );
    if (atom_types == 2)
      std::fprintf(file, "2\t%.8f\n", mass_);

    // Atoms, id type x y z
    std::fprintf(file, "\nAtoms\n\n");
    for (int i=0; i<node_positions_.size(); i++)
      std::fprintf(file, "%i\t%i\t%.8f\t%.8f\t%.8f\n", i, get_type(node_if_emphasis_[i]),
          node_positions_[i][0], node_positions_[i][1], node_positions_[i][2]);

    // Bonds, id type a b
    std::fprintf(file, "\nBonds\n\n");
    int i = 0;
    for (; i<bond_relations1_.size(); i++)
      std::fprintf(file, "%i\t1\t%i\t%i\n", i, bond_relations1_[i][0], bond_relations1_[i][1]);
    for (; i<bond_relations1_.size()+bond_relations2_.size(); i++)
      std::fprintf(file, "%i\t2\t%i\t%i\n", i, bond_relations2_[i][0], bond_relations2_[i][1]);
  }
  std::fclose(file);

  // dump trajectory, append file when time step is not 0
  file = std::fopen((file_name + ".dump").c_str(), (__time_step != 0) ? "a" : "w");
    
  // header
  std::fprintf(file, __dump_file_header, 
    __time_step, node_positions_.size(), 0., boundary_max, 0., boundary_max, 0., boundary_max
  );
  if (DUMP_CHECK(Metadata::kPrintVelocities, dump_type))
    std::fprintf(file, " vx vy vz");
  
  // self-defined contents, in kDumpMetaData
  for (auto i : Metadata::kDumpMetaData) {
    if (DUMP_CHECK(i.dump_type, dump_type))
      std::fprintf(file, " %s", i.name);
  }
  std::fprintf(file, "\n");

  // data
  for (int i=0; i<node_positions_.size(); i++) {
    // basic: id type xs ys zs
    std::fprintf(file, "%i\t%i\t%.8f\t%.8f\t%.8f", i, get_type(node_if_emphasis_[i]),
      node_positions_[i][0], node_positions_[i][1], node_positions_[i][2]);
    if (DUMP_CHECK(Metadata::kPrintVelocities, dump_type))
      std::fprintf(file, "\t%.8f\t%.8f\t%.8f", node_velocities_[i][0], node_velocities_[i][1],
        node_velocities_[i][2]);
    
    auto a = h_get_positions(i, node_adjacents_bonds1_[i]);
    auto b = h_get_positions(i, node_adjacents_bonds2_[i]);
    auto c = h_get_positions(i, node_adjacents_curvature_[i]);

    for (auto j : Metadata::kDumpMetaData)
      if (DUMP_CHECK(j.dump_type, dump_type))
        std::fprintf(file, "\t%.8f", j.func_num(*this, a, b, c));
    std::fprintf(file, "\n");
  }
  std::fclose(file);

  return;
}

void ModelSystem::store(std::string file_name) {

}

void ModelSystem::load(std::string file_name) {

}

void ModelSystem::update(bool just_velocity) {
  // gradients of curvature
  Kokkos::View<CoreMath::Array<CoreMath::Vector>*> gradients("", node_velocities_.size());

  // update using bonds
  Kokkos::parallel_for(node_velocities_.size(), KOKKOS_CLASS_LAMBDA(const int i) {
    // if it's not boundary or next to bound, curvature will take effect
    if (!node_if_rigid1_(i) && !node_if_rigid2_(i) && 
        !node_if_next_to_rigid1_(i) && !node_if_next_to_rigid2_(i)) {
      auto positions = d_get_positions(i, node_adjacents_curvature_(i));
      gradients(i) = curvature_gradient(positions);
    } else {
      gradients(i) = CoreMath::Array<CoreMath::Vector>(node_adjacents_curvature_(i).size());
    }

    CoreMath::Vector reduced;
    CoreMath::Array<int> position1 = node_adjacents_bonds1_(i), 
        position2 = node_adjacents_bonds2_(i);
    // If it's boundary nodes, don't count its related boundary nodes (for rigid body).
    if (node_if_rigid1_(i) || node_if_rigid2_(i)) {
      for (auto j=position1.begin(); j!=position1.end(); ) {
        // The node has been erased, the same place has the diffrent node needing a second count.
        if (node_if_rigid1_(*j) || node_if_rigid2_(*j))
          position1.erase(j);
        else j++;
      }
      for (auto j=position2.begin(); j!=position2.end(); ) {
        if (node_if_rigid1_(*j) || node_if_rigid2_(*j))
          position2.erase(j);
        else j++;
      }
    }
    // total force arised from bonds of type 1
    for (auto j : d_get_positions(i, position1))
      reduced += bond1_gradient(j);
    // total force arised from bonds of type 2
    for (auto j : d_get_positions(i, position2))
      reduced += bond2_gradient(j);

    // here node velocities are just -div(), divide by damp_coeff_ later
    node_velocities_(i) = reduced;
  });

  // another loop because we need to wait for every gradient finish
  Kokkos::parallel_for(node_velocities_.size(), KOKKOS_CLASS_LAMBDA(const int i) {
    // force arised from other node's curvature
    CoreMath::Vector reduced;
    // reduced vector for this node itself, no need to use parallel_reduce
    for (int j=0; j<gradients(i).size(); j++)
      reduced += gradients(i)[j];

    for (int j=0; j<node_adjacents_curvature_(i).size(); j++) {
      int adj = node_adjacents_curvature_(i)[j];
      // find related bond of adj and i
      for (int k=0; k<node_adjacents_curvature_(adj).size(); k++)
        if (node_adjacents_curvature_(adj)[k] == i) {
          reduced += -gradients(adj)[k];
          break;
        }
    }

    // random number, avoid waste when temperature equals 0
    if (temperature_ != 0 && !node_if_rigid1_(i) && !node_if_rigid2_(i))
      reduced += rand_pool_.gen_vector(Kokkos::sqrt(2*damp_coeff_*temperature_*K_B/mass_));
    node_velocities_(i) = (node_velocities_(i) + reduced) / damp_coeff_;
  });

  node_velocities_.modify<MemorySpace>();
  Kokkos::fence();
  if (just_velocity)
    return;
  
  // Center of mass, inertia tensor, total force, total moment, angular acceleration
  CoreMath::Vector center1, tensor1, force1, moment1, center2, tensor2, force2, moment2;
  // update non-boundary and boundary nodes
  Kokkos::parallel_reduce(node_positions_.size(), KOKKOS_CLASS_LAMBDA(const int i,
      CoreMath::Vector& center_inner1, CoreMath::Vector& center_inner2) {
    if (node_if_rigid1_(i)) {
      center_inner1 += node_positions_(i);
    } else if (node_if_rigid2_(i)) {
      center_inner2 += node_positions_(i);
    } else {
      node_positions_(i) += node_velocities_(i) * step_length_;
    }
  }, center1, center2);

  center1 = center1 / node_if_rigid1_count_;
  center2 = center2 / node_if_rigid2_count_;

  // calculate with rigid body
  Kokkos::parallel_reduce(node_positions_.size(), KOKKOS_CLASS_LAMBDA(const int i,
      CoreMath::Vector& force_inner1, CoreMath::Vector& force_inner2,
      CoreMath::Vector& moment_inner1, CoreMath::Vector& moment_inner2,
      CoreMath::Vector& tensor_inner1, CoreMath::Vector& tensor_inner2) {
    if (node_if_rigid1_(i)) {
      auto force = damp_coeff_ * node_velocities_(i);
      force_inner1 += force;
      auto t = node_positions_(i) - center1;
      moment_inner1 += CoreMath::cross(t, force);
      tensor_inner1 += CoreMath::Vector(t[1]*t[1]+t[2]*t[2], t[0]*t[0]+t[2]*t[2], t[0]*t[0]+t[1]*t[1]);
    } else if (node_if_rigid2_(i)) {
      auto force = damp_coeff_ * node_velocities_(i);
      force_inner2 += force;
      auto t = node_positions_(i) - center2;
      moment_inner2 += CoreMath::cross(t, force);
      tensor_inner2 += CoreMath::Vector(t[1]*t[1]+t[2]*t[2], t[0]*t[0]+t[2]*t[2], t[0]*t[0]+t[1]*t[1]);
    }
  }, force1, force2, moment1, moment2, tensor1, tensor2);

  // principal axis approximation, we don't need precise handle for boundary
  tensor1 = node_if_rigid1_count_ * CoreMath::Vector(moment1[0]/tensor1[0], 
      moment1[1]/tensor1[1], moment1[2]/tensor1[2]);
  tensor2 = node_if_rigid2_count_ * CoreMath::Vector(moment2[0]/tensor2[0], 
      moment2[1]/tensor2[1], moment2[2]/tensor2[2]);
  
  Kokkos::parallel_for(node_positions_.size(), KOKKOS_CLASS_LAMBDA(const int i) {
    if (node_if_rigid1_(i)) {
      auto t = node_positions_(i) - center1;
      node_positions_(i) += (force1 + CoreMath::cross(tensor1, t)) *
          step_length_ / damp_coeff_;
    } else if (node_if_rigid2_(i)) {
      auto t = node_positions_(i) - center2;
      node_positions_(i) += (force2 + CoreMath::cross(tensor2, t)) *
          step_length_ / damp_coeff_;
    }
  });
  
  __time_step++;
  node_positions_.modify<MemorySpace>();
}
