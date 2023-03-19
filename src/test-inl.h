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

#include <Kokkos_Core.hpp>

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

/**
 * @test test the result of parallel
 */
void test_parallel() {
  Kokkos::initialize(); {
  using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
  using HostMirrorSpace = Kokkos::Impl::HostMirror<MemorySpace>::Space;

  auto test_adj = test_adjacents();
  double para[3] = {1, Kokkos::sqrt(3)/2, 0};

  // init a list of nodes for test
  CoreMath::View<CoreMath::Array<CoreMath::Vector>> adjacents;
  CoreMath::View<CoreMath::Vector> results;
  adjacents.init(1000000);
  results.init(adjacents.size());
  for (int i=0; i<adjacents.size(); i++)
    adjacents[i] = test_adj;

  adjacents.modify<HostMirrorSpace>();
  adjacents.sync<MemorySpace>();

  auto t1 = std::chrono::high_resolution_clock::now();
  // non-parallel compute
  CoreMath::Vector* result = new CoreMath::Vector[adjacents.size()];
  for (int i=0; i<adjacents.size(); i++) {
    result[i] = CoreMath::Vector();
    auto inner = CoreEnergy::curvature_gradient(para, test_adj);
    for (auto j : inner)
      result[i] += j;
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  // test for parallel
  Kokkos::parallel_for(adjacents.size(), KOKKOS_LAMBDA(const int i) {
    CoreMath::Array<CoreMath::Vector> inner = 
      CoreEnergy::curvature_gradient(para, adjacents(i));

    results(i) = CoreMath::Vector();
    for (auto j : inner)
      results(i) += j;
  });

  Kokkos::fence(); // avoid asynchronous
  auto t3 = std::chrono::high_resolution_clock::now();

  results.modify<MemorySpace>();
  results.sync<HostMirrorSpace>();
  Kokkos::fence(); // avoid asynchronous
  auto t4 = std::chrono::high_resolution_clock::now();

  // count how much time is used
  std::chrono::duration<double> duration1 = t2 - t1, duration2 = t3 - t2, duration3 = t4 - t3;
  std::printf("total time: \tnon-parallel=%fms \tparallel=%fms \tmemory copy=%fms\n", 
      duration1.count()*1000, duration2.count()*1000, duration3.count()*1000);
  
  std::printf("given result: \tx=%.16f \ty=%.16f \tz=%.16f\n", 
    results[0][0], results[0][1], results[0][2]);
  std::printf("actual result: \tx=%.16f \ty=%.16f \tz=%.16f\n",
    result[0][0], result[0][1], result[0][2]);
  
  for (int i=0; i<adjacents.size(); i++) {
    if ((results[i] - result[i])[0] > 1e-12 || (results[i] - result[i])[0] < -1e-12) {
      std::printf("error: thread %i failed\n", i);
      break;
    }
  }
  delete []result;
  } Kokkos::finalize();
}