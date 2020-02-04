cmake_minimum_required(VERSION 3.13)
project(exp)
project(main)

set(CMAKE_CXX_STANDARD 14)

add_executable(exp matrix_utils.h runc_m.h exp.cpp)
add_executable(main matrix_utils.h main.cpp)
