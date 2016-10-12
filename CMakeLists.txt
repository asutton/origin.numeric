# Copyright (c) 2008-2016 Andrew Sutton
# All rights reserved

cmake_minimum_required(VERSION 3.0)
project(origin.numeric CXX)

enable_testing()

set(CMAKE_CXX_FLAGS "-std=c++1z -fconcepts")

include_directories(.)

# TODO: Move this into a position where it could be reused.
macro(add_unit_test target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} graph)
  add_test(${target} ${target})
endmacro()

add_subdirectory(origin)
