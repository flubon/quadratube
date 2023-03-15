/**
 * @file test-inl.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief test
 * @version 0.0.1
 * @date 2023-03-14
 * 
 * @copyright Copyright (c) 2023
 */

#include <stdio.h>

#include <chrono>

#include "core/math.h"
#include "core/energy.h"

/// @brief example adjacents for check
/// @return
inline CoreMath::Array<CoreMath::Vector> test_adjacents() {
  CoreMath::Array<CoreMath::Vector> example(6);
  example[0] = CoreMath::Vector(-0.2848985, -0.54609, 0.7857);
  example[1] = CoreMath::Vector(-0.3673606, -0.90738, -0.1429);
  example[2] = CoreMath::Vector(-0.192107, -0.31690, -0.9286);
  example[3] = CoreMath::Vector(0.442461, 0.42851, -0.7857);
  example[4] = CoreMath::Vector(0.765357, 0.61036, 0.1428);
  example[5] = CoreMath::Vector(0.249141, 0.27434, 0.9285);
  return example;
}

/**
 * @test compute whether the numerical derivation result is the same 
 *     with the analytic derivation result
 */
void test_derivative() {
  double test_bending_rigidity = 10, precision = 1e-9;
  CoreMath::Vector numerical_result, analytic_result;

  auto others = test_adjacents();

  auto t1 = std::chrono::high_resolution_clock::now();
  // begin numerical derivation
  for (int i=0; i<1000; i++) {
    auto other1(others), other2(others), other3(others);
    for (int i=0; i<others.size(); i++) {
      other1[i] += CoreMath::Vector(precision, 0, 0);
      other2[i] += CoreMath::Vector(0, precision, 0);
      other3[i] += CoreMath::Vector(0, 0, precision);
    }
    double origin = CoreEnergy::curvature_energy(&test_bending_rigidity, others);
    numerical_result = CoreMath::Vector(
      origin - CoreEnergy::curvature_energy(&test_bending_rigidity, other1),
      origin - CoreEnergy::curvature_energy(&test_bending_rigidity, other2),
      origin - CoreEnergy::curvature_energy(&test_bending_rigidity, other3)
    ) / precision;
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  // begin analytic derivation
  for (int i=0; i<1000; i++) {
    auto part_result = CoreEnergy::curvature_gradient(&test_bending_rigidity, others);
    analytic_result = CoreMath::Vector();
    for (auto i : part_result)
      analytic_result += -i;
  }
  
  auto t3 = std::chrono::high_resolution_clock::now();
  // count how much time is used
  std::chrono::duration<double> duration1 = t2 - t1, duration2 = t3 - t2;
  
  std::printf("numerical result: \tx=%.16f \ty=%.16f \tz=%.16f\n", 
      numerical_result[0], numerical_result[1], numerical_result[2]);
  std::printf("analytic result: \tx=%.16f \ty=%.16f \tz=%.16f\n",
      analytic_result[0], analytic_result[1], analytic_result[2]);
  std::printf("total time: \tnumerical=%fms \tanalytic=%fms\n", 
      duration1.count()*1000, duration2.count()*1000);
}