/**
 * @file initializer3.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "model/initializer.h"

#include <Kokkos_Core.hpp>

namespace ModelInitializer {

int Initializer::dual(ObjectType tp, int ck, int a, int b) {
  int il = -1, ir = -1;
  for (auto i : __system.node_adjacents_bonds1_[b]) {
    auto j = __system.node_adjacents_bonds1_[a].find(i);
    if (j != __system.node_adjacents_bonds1_[a].end())
      ((il == -1) ? il : ir) = *j;
  }
  return (il == ck) ? ir : il;
}

void Initializer::remove(int node) {
  __system.node_positions_[node] = __system.node_positions_[__system.node_positions_.size()-1];
  __system.node_positions_.pop_back();
  // remove related bonds
  for (int i = __system.node_adjacents_bonds1_[node].size()-1; i>=0; i--)
    remove(Bond1, CoreMath::Pair<int>(node, __system.node_adjacents_bonds1_[node][i]));
  __system.node_adjacents_bonds1_[node] = __system.node_adjacents_bonds1_[__system.node_positions_.size()];
  __system.node_adjacents_bonds1_.pop_back();
  // keep consistency
  for (int i=0; i<__system.bond_relations1_.size(); i++)
    __system.bond_relations1_[i].replace(__system.node_positions_.size(), node);
  for (int i=0; i<__system.node_adjacents_bonds1_[node].size(); i++)
    *__system.node_adjacents_bonds1_[__system.node_adjacents_bonds1_[node][i]]
      .find(__system.node_positions_.size()) = node;
}

/**
 * @brief 
 * @details Rember we need to init such objects:
 * 1. node_positions_: will be used many times
 * 2. node_velocities_: calculate by update
 * 3-4. node_adjacents_bonds1_ & node_adjacents_curvature_: 
 *   keep the same in memory space (won't be used by host)
 * 5. node_adjacents_bonds2_: just need to initialize on memory space
 * 6-7. node_if_emphasis_ & node_if_emphasis_count_: use begin and end
 * 8-11. node_if_rigid1_ & node_if_next_to_rigid1_ & node_if_rigid2_ & 
 *   node_if_next_to_rigid2_: check at the end
 * 12-15. node_if_rigid1_count_ & node_if_next_to_rigid1_count_ & node_if_rigid2_count_
 *   & node_if_next_to_rigid2_count_: calculate by structures before
 * 16. bond_relations1_: used by host, don't need to sync
 * 17. bond_relations2_: won't be used at all, don't need to care
 * 
 * To seeking security, we use `INIT 1`, `CHANGE 1`, `FINISH 1` to mark status.
 * 
 * @param init_para
 */
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 1. generate perfect model's positions
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  double m = static_cast<double>(init_para.m),
      n = static_cast<double>(init_para.n);

  // angle A, B (between x axis and m, n), radius (normalized)
  double A = (2*n == m) ? PI / 2 : Kokkos::atan(Kokkos::sqrt(3)/2 * m / (n - m/2));
  double B = (2*m == n) ? PI / 2 : Kokkos::atan(Kokkos::sqrt(3)/2 * n / (m - n/2));
  double r = 1/PI/2 * Kokkos::sqrt(m*m+n*n-m*n);

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
      auto near_then_push_all = [=](int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j)][2] -
            __system.node_positions_[a][2]) < 2*init_para.rest_len) {
          __system.bond_relations1_.push_back({flat(i, j), a});
          __system.node_adjacents_bonds1_[flat(i, j)].push_back(a);
        }
      };
      // don't push_back bond again because it has been counted
      auto near_then_push = [=](int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j)][2] -
            __system.node_positions_[a][2]) < 2*init_para.rest_len)
          __system.node_adjacents_bonds1_[flat(i, j)].push_back(a);
      };
      // CHANGE 3. CHANGE 16.
      near_then_push_all(flat(i+1, j));
      near_then_push_all(flat(i, j+1));
      near_then_push_all(flat(i-1, j+1));
      near_then_push(flat(i-1, j));
      near_then_push(flat(i, j-1));
      near_then_push(flat(i+1, j-1));
    }
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 3. init checkpoints according to direction etc.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  // dislocations 1-4: begin[0], begin[1], end[0], end[1]
  int dislocations[4] = {0};
  dislocations[1] = dislocations[2] = flat(0, (init_para.repeat+1)*init_para.m / 2) + init_para.bn;
  CoreMath::Array<int>& center = __system.node_adjacents_bonds1_[dislocations[2]];
  if (init_para.direction > 0)
    dislocations[0] = dislocations[3] = center[Kokkos::abs(init_para.direction)];
  if (init_para.direction < 0)
    dislocations[0] = dislocations[3] = 
        center[init_para.direction == -1 ? 5 : (Kokkos::abs(init_para.direction) - 2)];
  int checkpoint = rot(Bond1, center[Kokkos::abs(init_para.direction) - 1], 
      dislocations[2], dislocations[3]);
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 4. add glides
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  for (int i=0; i < init_para.glide; i++) {
    // quadrilateral of glide
    int new_end0 = rot(Bond1, checkpoint, dislocations[2], dislocations[3]);
    int new_end1 = rot(Bond1, dislocations[2], new_end0, dislocations[3]);

    remove(Bond1, CoreMath::Pair<int>(dislocations[3], new_end0));
    insert(Bond1, CoreMath::Pair<int>(dislocations[2], new_end1), 
        between(Bond1, dislocations[2], dislocations[3], new_end0), 
        between(Bond1, new_end1, dislocations[3], new_end0));

    // recover
    checkpoint = dislocations[2];
    dislocations[2] = new_end0;
    dislocations[3] = new_end1;
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 5. add climbs with nodes decreased
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  if (init_para.climb < 0)
  for (int i=0; i<-init_para.climb; i++) {
    // hexagon of climb
    int center = rot(Bond1, checkpoint, dislocations[2], dislocations[3]);
    int c1 = rot(Bond1, dislocations[2], center, dislocations[3]);
    int c2 = rot(Bond1, dislocations[3], center, c1);
    int c3 = rot(Bond1, c1, center, c2);
    int new_end0 = rot(Bond1, c2, center, c3);
        
    remove(center);
    insert(Bond1, CoreMath::Pair<int>(dislocations[2], c1), 
        between(Bond1, dislocations[2], dislocations[3], new_end0), 
        between(Bond1, c1, dislocations[3], c2));
    insert(Bond1, CoreMath::Pair<int>(dislocations[2], c2), 
        between(Bond1, dislocations[2], c1, new_end0), 
        between(Bond1, c2, c1, c3));
    insert(Bond1, CoreMath::Pair<int>(dislocations[2], c3), 
        between(Bond1, dislocations[2], c2, new_end0), 
        between(Bond1, c3, c2, new_end0));
        
    dislocations[3] = dislocations[2];
    dislocations[2] = new_end0;
    checkpoint = rot(Bond1, c3, dislocations[2], dislocations[3]);
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 6. add climbs with nodes increased
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  if (init_para.climb > 0)
  for (int i=0; i<init_para.climb; i++) {
    // pentagon of climb
    int c1 = rot(Bond1, checkpoint, dislocations[2], dislocations[3]);
    int c2 = rot(Bond1, dislocations[2], c1, dislocations[3]);
    int new_end1 = rot(Bond1, c1, c2, dislocations[3]);
        
    remove(Bond1, CoreMath::Pair<int>(dislocations[3], c1));
    remove(Bond1, CoreMath::Pair<int>(dislocations[3], c2));
    int new_end0 = insert((__system.node_positions_[dislocations[3]] + 
        __system.node_positions_[c2] + __system.node_positions_[c1])/3);
    // don't need between because it's order is fixed
    insert(Bond1, CoreMath::Pair<int>(new_end0, dislocations[2]),
        __system.node_adjacents_bonds1_[new_end0].end(),
        between(Bond1, dislocations[2], dislocations[3], c1));
    insert(Bond1, CoreMath::Pair<int>(new_end0, dislocations[3]),
        __system.node_adjacents_bonds1_[new_end0].end(), 
        between(Bond1, dislocations[3], dislocations[2], new_end1));
    insert(Bond1, CoreMath::Pair<int>(new_end0, new_end1), 
        __system.node_adjacents_bonds1_[new_end0].end(),
        between(Bond1, new_end1, dislocations[3], c2));
    insert(Bond1, CoreMath::Pair<int>(new_end0, c2),
        __system.node_adjacents_bonds1_[new_end0].end(),
        between(Bond1, c2, new_end1, c1));
    insert(Bond1, CoreMath::Pair<int>(new_end0, c1),
        __system.node_adjacents_bonds1_[new_end0].end(),
        between(Bond1, c1, dislocations[2], c2));
        
    checkpoint = dislocations[3];
    dislocations[2] = new_end0;
    dislocations[3] = new_end1;
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 7. set up emphasis for dump and sync
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  __system.node_positions_.modify_host();
  __system.node_positions_.sync_device();
  __system.node_adjacents_bonds1_.modify_host();
  __system.node_adjacents_bonds1_.sync_device();

  __system.node_if_emphasis_.init(__system.node_positions_.size(), nodes_number);
  __system.node_if_emphasis_count_ = 4;
  for (int i=0; i<__system.node_positions_.size(); i++)
    if (i == dislocations[0] || i == dislocations[1] || 
        i == dislocations[2] || i == dislocations[3])
      __system.node_if_emphasis_[i] = true;
  
  __system.node_if_emphasis_.modify_host();
  __system.node_if_emphasis_.sync_device();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 8. set up boundary nodes for update
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  __system.node_if_rigid1_.init(__system.node_positions_.size());
  __system.node_if_rigid2_.init(__system.node_positions_.size());

  // FINISH 8. FINISH 10. FINISH 12. FINISH 14. boundary count
  Kokkos::parallel_reduce(__system.node_positions_.size(), 
    KOKKOS_CLASS_LAMBDA(const int i, int& inner1, int& inner2){
    // if it's a boundary node
    if (__system.node_adjacents_bonds1_(i).size() != 6 && !__system.node_if_emphasis_(i)) {
      // if it's next to bottom
      if (__system.node_positions_(i)[2] < 5*init_para.rest_len) {
        __system.node_if_rigid1_(i) = true;
        inner1++;
        __system.node_if_rigid2_(i) = false;
      } else {
        __system.node_if_rigid2_(i) = true;
        inner2++;
        __system.node_if_rigid1_(i) = false;
      }
    } else {
      __system.node_if_rigid1_(i) = false;
      __system.node_if_rigid2_(i) = false;
    }
  }, __system.node_if_rigid1_count_, __system.node_if_rigid2_count_);

  // for consistent with function store()
  __system.node_if_rigid1_.modify_device();
  __system.node_if_rigid1_.sync_host();
  __system.node_if_rigid2_.modify_device();
  __system.node_if_rigid2_.sync_host();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 9. set up nodes for boundary for update
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  __system.node_if_next_to_rigid1_.init(__system.node_positions_.size());
  __system.node_if_next_to_rigid2_.init(__system.node_positions_.size());

  // FINISH 9. FINISH 11. FINISH 13. FNISH 15. check next to boundary
  Kokkos::parallel_reduce(__system.node_positions_.size(), 
    KOKKOS_CLASS_LAMBDA(const int i, int& inner1, int& inner2){
    // if it's not a boundary node
    if (!__system.node_if_rigid1_(i) && !__system.node_if_rigid2_(i)) {
      for (auto j : __system.node_adjacents_bonds1_(i)) {
        // if it's next to boundary node
        if (__system.node_if_rigid1_(i)) {
          __system.node_if_next_to_rigid1_(i) = true;
          inner1++;
          __system.node_if_next_to_rigid2_(i) = false;
          // skip additional assignment (equal to break; and continue; immediately)
          goto ContinueLoop;
        } else if (__system.node_if_rigid2_(i)) {
          __system.node_if_next_to_rigid2_(i) = true;
          inner2++;
          __system.node_if_next_to_rigid1_(i) = false;
          goto ContinueLoop;
        }
      }
      // if it's not next to boundary node
      __system.node_if_next_to_rigid1_(i) = false;
      __system.node_if_next_to_rigid2_(i) = false;
      ContinueLoop: ;
    }
  }, __system.node_if_next_to_rigid1_count_, __system.node_if_next_to_rigid2_count_);

  __system.node_if_next_to_rigid1_.modify_device();
  __system.node_if_next_to_rigid1_.sync_host();
  __system.node_if_next_to_rigid2_.modify_device();
  __system.node_if_next_to_rigid2_.sync_host();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//   STEP 10. set redundant objects for model3
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  // FINISH 5. `node_adjacents_bonds2_` wil be used by update, so make sure it has been initialized.
  __system.node_adjacents_bonds2_.init(__system.node_positions_.size());
  for (int i=0; i<__system.node_positions_.size(); i++)
    __system.node_adjacents_bonds2_[i] = CoreMath::Array<int>();
  __system.node_adjacents_bonds2_.modify_host();
  __system.node_adjacents_bonds2_.sync_device();

  // FINISH 4. `node_adjacents_curvature_` must keep consistent with `node_adjacents_bonds1_`, but 
  //    won't be used by host. So just copy MemorySpace.
  __system.node_adjacents_curvature_.init(__system.node_positions_.size());
  Kokkos::deep_copy(__system.node_adjacents_curvature_.view_device(),
      __system.node_adjacents_bonds1_.view_device());
  __system.node_adjacents_bonds1_.modify_device();
  __system.node_adjacents_bonds1_.sync_host();
  // FINISH 17. `bond_relations2_` won't be used at all.

  // FINISH 2. set velocities
  __system.node_velocities_.init(__system.node_positions_.size());
  __system.update(true);
}

} // namespace ModelInitializer