cmake_minimum_required(VERSION 3.13)
project(comp)
project(main)
project(dev)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS "-pg -Wall -O3")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(main matrix_utils.h AO2_trajectory.h AO2.h main.cpp)
add_executable(dev matrix_utils.h AO2_trajectory.h dev.cpp)
