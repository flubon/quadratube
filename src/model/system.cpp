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
  sync<HostMirrorSpace>();

  // boundary: x, y, z: [0, 30*bond1_rest_length_]
  double boundary_max = 30*bond1_rest_length_;

  int atom_types = DUMP_CHECK(Metadata::kPrintDislocations, dump_type) ? 2 : 1;
  // default type is 1, dislocations use type 2
  int* types = new int[node_positions_.size()];
  for (auto i=0; i<node_positions_.size(); i++)
    types[i] = 1;
  if (atom_types == 2)
    for (int i=0; i<part_node_emphasis_.size(); i++)
      types[part_node_emphasis_[i]] = 2;

  // print bonds into data file, skip if file exists already.
  FILE *file = std::fopen((file_name + ".data").c_str(), "r");
  if (DUMP_CHECK(Metadata::kPrintBond, dump_type) || file == NULL) {
    file = std::fopen((file_name + ".data").c_str(), "w");

    // header
    size_t bonds = DUMP_CHECK(Metadata::kExcludeBondType2, dump_type) ?
        bond_relations1_.size() : bond_relations1_.size() + bond_relations2_.size();
    int bond_types = DUMP_CHECK(Metadata::kExcludeBondType2, dump_type) ? 1 : 2;
    std::fprintf(file, __data_file_header,
      node_positions_.size(), bonds, atom_types, bond_types,          // atom & bond number
      0., boundary_max, 0., boundary_max, 0., boundary_max, this->mass_  // boundary & masses
    );
    if (atom_types == 2)
      std::fprintf(file, "2\t%.8f\n", this->mass_);

    // Atoms, id type x y z
    std::fprintf(file, "\nAtoms\n\n");
    for (int i=0; i<node_positions_.size(); i++)
      std::fprintf(file, "%i\t%i\t%.8f\t%.8f\t%.8f\n", i, types[i],
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
    std::fprintf(file, "%i\t%i\t%.8f\t%.8f\t%.8f", i, types[i],
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

  delete[] types;
  return;
}

void ModelSystem::store(std::string file_name) {

}

void ModelSystem::load(std::string file_name) {

}

void ModelSystem::update() {
  CoreMath::View<CoreMath::Array<CoreMath::Vector>> curvature_gradients;

}
