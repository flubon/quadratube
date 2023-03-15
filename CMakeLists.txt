cmake_minimum_required(VERSION 3.16)
project(quadratube VERSION 0.0.1)

# 文件头
include_directories(src)

# Kokkos 文件夹，将第一个替换成自己的文件
set(Kokkos_DIR /home/bovera/libkokkos/lib/cmake/Kokkos)
find_package(Kokkos REQUIRED)

# 设置输出为调试，不需要自行注释掉
set(CMAKE_BUILD_TYPE "Debug")

# 源文件，编译和链接成可执行文件
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(quadratube ${SOURCES})
target_link_libraries(quadratube Kokkos::kokkos)
