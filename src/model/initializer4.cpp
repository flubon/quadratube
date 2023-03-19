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

#if MODEL_TYPE == 4

namespace ModelInitializer {

void Initializer4::init(Parameters4 init_para) {
  __system.init_all(1000, 3000);
}

} // namespace ModelInitializer

#endif // #if MODEL_TYPE == 3