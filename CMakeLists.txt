cmake_minimum_required(VERSION 3.16)
project(quadratube VERSION 0.0.1)

# triangular or quadrangular
option(Model3 "triangular model" OFF)
option(DebugType "cmake debug build type" ON)

# header files' root path
include_directories(src)

# Kokkos directory, replace to your path
set(Kokkos_DIR /home/bovera/libkokkos/lib/cmake/Kokkos)
find_package(Kokkos REQUIRED)

# Debug mode, no optimization, can use gdb
if(DebugType)
  set(CMAKE_BUILD_TYPE "Debug")
else()
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Source files, compile and link into executable
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Choose corresponding file and macro to use
if(Model3)
  list(REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/model/initializer4.cpp)
  add_definitions(-DMODEL_TYPE=3)
else()
  list(REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/model/initializer3.cpp)
  add_definitions(-DMODEL_TYPE=4)
endif()

add_executable(quadratube ${SOURCES})
target_link_libraries(quadratube Kokkos::kokkos)

# Force to use c++17 and avoid incompatibility between compilers
set_target_properties(quadratube 
  PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)