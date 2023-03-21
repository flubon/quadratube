/**
 * @file initializer4.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "model/initializer.h"

namespace ModelInitializer {

int Initializer::dual(ObjectType tp, int ck, int a, int b) {
  int il = -1, ir = -1;
  for (auto i : node_adjacents(tp)[b]) {
    auto j = node_adjacents(tp)[a].find(i);
    if (j != node_adjacents(tp)[a].end())
      ((il == -1) ? il : ir) = *j;
  }
  return (il == ck) ? ir : il;
}

void Initializer::remove(int node) {}

void Initializer::init(Parameters init_para) {
  // this parameters struct is used by other 'small' functions
  __para = init_para;
  // allocate size of system
  int nodes_number = init_para.m * init_para.n * init_para.repeat +
      ((init_para.climb < 0) ? 0 : init_para.climb);
  // INIT 1. INIT 3. INIT 16.
  __system.node_positions_.init(0, nodes_number);
  // node numbers of perfect model, avoid resizing after perfect model construct
  __system.node_adjacents_bonds1_.init(init_para.m*init_para.n*init_para.repeat, nodes_number);
  __system.bond_relations1_.init(0, 3 * nodes_number);
  __system.node_adjacents_bonds2_.init(init_para.m*init_para.n*init_para.repeat, nodes_number);
  __system.bond_relations2_.init(0, 3 * nodes_number);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 1. generate perfect model's positions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  double m = static_cast<double>(init_para.m),
      n = static_cast<double>(init_para.n);

  // angle A, B (between x axis and m, n), radius (normalized)
  double A = (n == 0) ? PI : Kokkos::atan(m / n);
  double B = PI/2 - A;
  double r = Kokkos::sqrt(m*m+n*n)/PI/2;

  for (int j = 0; j < init_para.repeat*init_para.m; j++)
    for (int i = 0; i < init_para.n; i++) {
      double x = -(i-n)*Kokkos::cos(A) + j*Kokkos::cos(B);
      double z = (i-n)*Kokkos::sin(A) + j*Kokkos::sin(B);
      // if below the axis
      if (z < 0) {
        x += init_para.repeat*m*Kokkos::cos(B);
        z += init_para.repeat*m*Kokkos::sin(B);
      }

      // transfrom to 3-dimensional coordinates
      double y = r * (Kokkos::sin(x/r) + 1);
      x = r * (Kokkos::cos(x/r) + 1);
      // CHANGE 1. add to queue
      __system.node_positions_.push_back(CoreMath::Vector(init_para.rest_len*x,
        init_para.rest_len*y, init_para.rest_len*z));
    }
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 2. generate bonds and adjacents' relations of perfect model
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  for (int i = 0; i < init_para.n; i++)
    for (int j = 0; j < init_para.repeat*init_para.m; j++) {
      auto near_then_push_all = [=](ObjectType tp, int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j)][2] -
            __system.node_positions_[a][2]) < 2*init_para.rest_len) {
          bond_relations(tp).push_back({flat(i, j), a});
          node_adjacents(tp)[flat(i, j)].push_back(a);
        }
      };
      // don't push_back bond again because it has been counted
      auto near_then_push = [=](ObjectType tp, int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j)][2] -
            __system.node_positions_[a][2]) < 2*init_para.rest_len)
          node_adjacents(tp)[flat(i, j)].push_back(a);
      };
      // CHANGE 3. CHANGE 16.
      near_then_push_all(Bond1, flat(i+1, j));
      if ((i+j)%2 == 0)
          near_then_push_all(Bond2, flat(i+1, j+1));
      near_then_push_all(Bond1, flat(i, j+1));
      if ((i+j)%2 == 0)
        near_then_push_all(Bond2, flat(i-1, j+1));
      near_then_push(Bond1, flat(i-1, j));
      if ((i+j)%2 == 0)
        near_then_push(Bond2, flat(i-1, j-1));
      near_then_push(Bond1, flat(i, j-1));
      if ((i+j)%2 == 0)
        near_then_push(Bond2, flat(i+1, j-1));
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 3. init checkpoints according to direction etc.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  // dislocations 1-4: begin[0], begin[1], end[0], end[1]
  int dislocations[4] = {0};
  dislocations[1] = dislocations[2] = flat(0, (init_para.repeat+1)*init_para.m / 2) + init_para.bn;
  CoreMath::Array<int>& center = __system.node_adjacents_bonds1_[dislocations[2]];

  // TODO, two directions
  int il = center[0], ir = center[1];
  int checkpoint = center[0];
  dislocations[0] = dislocations[3] = dual(Bond1, dislocations[2], il, ir);
  insert(Bond1, CoreMath::Pair<int>(dislocations[2], dislocations[3]), 
      between(Bond1, dislocations[2], il, ir), 
      between(Bond1, dislocations[3], il, ir));
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 4. add glides
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  for (int i=0; i < init_para.glide; i++) {
    // quadrilateral of glide
    int new_end1 = rot(Bond1, checkpoint, dislocations[2], dislocations[3]);
    int new_end0 = rot(Bond1, dislocations[3], dislocations[2], new_end1);
    int c0 = dual(Bond1, dislocations[2], new_end1, new_end0);

    remove(Bond1, CoreMath::Pair<int>(dislocations[2], new_end1));
    insert(Bond1, CoreMath::Pair<int>(new_end0, new_end1), 
        between(Bond1, new_end0, dislocations[3], c0), 
        between(Bond1, new_end1, dislocations[2], c0));

    // recover
    checkpoint = dislocations[2];
    dislocations[2] = new_end0;
    dislocations[3] = new_end1;
  }

}

} // namespace ModelInitializer