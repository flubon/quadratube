/**
 * @file initializer.cpp
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

void Initializer3::init(Parameters3 init_parameter) {

}

void Initializer4::init(Parameters4 init_parameter) {
  __system.init(1000, 1000, 1000);
}

} // namespace ModelInitializer