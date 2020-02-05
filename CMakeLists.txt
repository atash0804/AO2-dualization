cmake_minimum_required(VERSION 3.13)
project(exp)
project(main)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS "-pg -Wall")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")


add_executable(exp matrix_utils.h runc_m.h exp.cpp)
add_executable(main matrix_utils.h main.cpp)
