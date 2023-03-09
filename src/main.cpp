/**
 * @file main.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief Main function
 * @version 0.0.1
 * @date 2023-03-07
 * 
 * @copyright Copyright (c) 2023
 */
#include <iostream>

#include "core/math.h"

int main() {
  CoreMath::Array<int> a(3);
  CoreMath::Array<int> b(a);
  b.push_back(4);
  std::cout << b[3];
  for (auto i : b) {
    std::cout << i << std::endl;
  }
}