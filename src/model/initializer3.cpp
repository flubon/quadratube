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

#if MODEL_TYPE == 3

namespace ModelInitializer {

namespace {

/// @brief calculate index from 2-dimensional coordinates
inline int flat(int i, int j, Parameters init_para) {
  j += init_para.m * (i/init_para.n);
  i %= init_para.n;
  if (i < 0) {
    i += init_para.n;
    j -= init_para.m;
  }
  j %= init_para.m * init_para.repeat;
  if (j < 0)
    j += init_para.m * init_para.repeat;
  return i + init_para.n*j;
}

} // namespace

void Initializer::init(Parameters init_para) {
  // allocate
  int nodes_number = init_para.m * init_para.n * init_para.repeat + 
      ((init_para.climb < 0) ? 0 : init_para.climb);
  __system.init_all(nodes_number, 3 * nodes_number);

  // begin perfect model generate
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
      // add to queue
      __system.node_positions_.push_back(CoreMath::Vector(init_para.rest_len*x,
        init_para.rest_len*y, init_para.rest_len*z));
    }
  
  __system.resize_nodes(__system.node_positions_.size());

  // generate bonds and adjacents' relations
  for (int i = 0; i < init_para.n; i++) 
    for (int j = 0; j < init_para.repeat*init_para.m; j++) {
      auto near_then_push_all = [=](int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j, init_para)][2] - __system.node_positions_[a][2]) 
            < 2*init_para.rest_len) {
          __system.bond_relations1_.push_back({flat(i, j, init_para), a});
          __system.node_adjacents_bonds1_[flat(i, j, init_para)].push_back(a);
        }
      };
      // don't push_back bond again because it has been counted
      auto near_then_push = [=](int a) {
        if (Kokkos::abs(__system.node_positions_[flat(i, j, init_para)][2] - __system.node_positions_[a][2])
            < 2*init_para.rest_len)
          __system.node_adjacents_bonds1_[flat(i, j, init_para)].push_back(a);
      };

      near_then_push_all(flat(i+1, j, init_para));
      near_then_push_all(flat(i, j+1, init_para));
      near_then_push_all(flat(i-1, j+1, init_para));
      near_then_push(flat(i-1, j, init_para));
      near_then_push(flat(i, j-1, init_para));
      near_then_push(flat(i+1, j-1, init_para));
    }
}

} // namespace ModelInitializer

#endif // if MODEL_TYPE == 3