/**
 * @file modifier.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-09
 * 
 * @copyright Copyright (c) 2023
 */
#ifndef QUADRATUBE_UTILS_MODIFIER_H_
#define QUADRATUBE_UTILS_MODIFIER_H_

#include "core/math.h"
#include "model/system.h"

namespace UtilsModifier {

/**
 * @brief modifier of system
 * 
 */
class Modifier {
  public:
    inline Modifier(ModelSystem& system): __system(system) {}

    /// @brief total energy in this box
    double total_energy(CoreMath::Vector range_l, CoreMath::Vector range_r);
    /// @brief total particle numbers in this box
    int total_particle(CoreMath::Vector range_l, CoreMath::Vector range_r);

  private:
    ModelSystem& __system;
};

} // namespace UtilsModifier

#endif // QUADRATUBE_UTILS_MODIFIER_H_