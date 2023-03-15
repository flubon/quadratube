/**
 * @file initializer.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-11
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef QUADRATUBE_MODEL_INITIALIZER_H_
#define QUADRATUBE_MODEL_INITIALIZER_H_

#include <string>

#include "core/math.h"
#include "model/system.h"

namespace ModelInitializer {

typedef struct {
  int m;  ///< parastichy number m
  int n;  ///< parastichy number n
  int repeat;  ///< 重复次数
  int direction; ///< 位错发生方向，正负号表示手性
  int glide;  ///< 滑移步数
  int climb;  ///< 攀移步数
  int bn; ///< begin position
} Parameters3;

typedef struct {
  int m;  //!< 叶序数 m
  int n;  //!< 叶序数 n
  int repeat;  //!< 重复次数
  int direction; //!< 位错发生方向，正负号表示手性
  int glide;  //!< 滑移步数
  int climb;  //!< 攀移步数
  int bn; //!< 位错在管上的初始位置
} Parameters4;

class Initializer3 {
  public:
    inline Initializer3(ModelSystem& system): __system(system) {}
    void init(Parameters3 init_parameter);

  private:
    ModelSystem& __system;
};

class Initializer4 {
  public:
    inline Initializer4(ModelSystem& system): __system(system) {}
    void init(Parameters4 init_parameter);

  private:
    ModelSystem& __system;
};

} // namespace ModelInitializer

#endif // QUADRATUBE_MODEL_INITIALIZER_H_
