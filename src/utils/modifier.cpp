/**
 * @file modifier.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-09
 * 
 * @copyright Copyright (c) 2023
 */
#include "utils/modifier.h"

namespace UtilsModifier {
  double Modifier::total_energy(CoreMath::Vector range_l, CoreMath::Vector range_r) {
    double count = 0;
    for (int i=0; i<__system.node_positions_.size(); i++) {
      auto others1 = __system.h_get_positions(i, __system.node_adjacents_bonds1_[i]);
      auto others2 = __system.h_get_positions(i, __system.node_adjacents_bonds2_[i]);
      auto others3 = __system.h_get_positions(i, __system.node_adjacents_curvature_[i]);
      for (int i=0; i<others1.size(); i++)
        count += __system.bond1_energy(others1[i]) / 2;
      for (int i=0; i<others2.size(); i++)
        count += __system.bond2_energy(others2[i]) / 2;
      count += __system.curvature_energy(others3);
    }
    return count;
  }

  int Modifier::total_particle(CoreMath::Vector range_l, CoreMath::Vector range_r) {
    int count = 0;
    for (int i=0; i<__system.node_positions_.size(); i++) {
      if (__system.node_positions_[i][0] > range_l[0] && __system.node_positions_[i][0] < range_r[0] &&
        __system.node_positions_[i][1] > range_l[1] && __system.node_positions_[i][1] < range_r[1] &&
        __system.node_positions_[i][2] > range_l[2] && __system.node_positions_[i][2] < range_r[2])
        count++;
    }
    return count;
  }
} // namespace UtilsModifier