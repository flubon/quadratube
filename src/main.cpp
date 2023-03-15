/**
 * @file main.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief Main function
 * @version 0.0.1
 * @date 2023-03-07
 * 
 * @copyright Copyright (c) 2023
 */

// #define RESTART

#include <stdio.h>

#include <Kokkos_Core.hpp>

#include "metadata.h"
#include "core/math.h"
#include "model/system.h"
#include "model/initializer.h"
#include "utils/modifier.h"

int main(int argc, char* argv[]) {
  // initialize kokkos in main function instead of class system
  // use {} to limit life cycle, avoiding `deallocate after Kokkos::finalize`
  Kokkos::initialize(argc, argv); {
  ModelSystem model;

  // set runtime parameters
  model.bond2_spring_constant_ = 1;
  model.curvature_bending_rigidity_ = 0.1;
  model.damp_coeff_ = 1;
  model.temperature_ = 0;

#ifdef RESTART
  model.load("restart.bin");
#else
  // initialize model with initializer
  ModelInitializer::Initializer4 initializer(model);
  ModelInitializer::Parameters4 parameters = {
    .m = 5, .n = 3, .repeat = 8, .direction = -1, .glide = 5, .climb = 0
  };
  initializer.init(parameters);
#endif

  // range of output box
  #define OUT_RANGE CoreMath::Vector(-INFINITY, -INFINITY, 29), CoreMath::Vector(INFINITY, INFINITY, 48)

  for (int k=10000; k>0; k--) {
    if (k%1000 == 0) {
      // dump current states
      model.dump("test", Metadata::kPrintAll);
      // use modifier to calculate related global quantities
      UtilsModifier::Modifier modifier(model);
      std::printf("total particle: %i, total energy: %.8f\n",
          modifier.total_particle(OUT_RANGE), modifier.total_energy(OUT_RANGE));
    }
    model.update();
  }

  model.store("restart.bin");
  } Kokkos::finalize(); // deconstruct before finalize

  return 0;
}