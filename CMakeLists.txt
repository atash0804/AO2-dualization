cmake_minimum_required(VERSION 3.13)
project(comp)
project(main)
project(dev)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS "-pg -Wall")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")


add_executable(comp matrix_utils.h runc_m.h comp.cpp)
add_executable(main matrix_utils.h main.cpp)
add_executable(dev matrix_utils.h dev.cpp)
