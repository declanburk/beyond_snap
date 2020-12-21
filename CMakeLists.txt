cmake_minimum_required(VERSION 3.0)

project(fast_snap)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DEIGEN_NO_DEBUG -std=c++11")

find_package(Eigen3 3.3 REQUIRED)
include_directories(include
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(
    test_fast_snap
    src/test_minimum_snap.cpp
)