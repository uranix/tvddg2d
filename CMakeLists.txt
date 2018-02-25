cmake_minimum_required(VERSION 3.0)

project(tvddg2d)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_CXX_FLAGS "-Wall --max-errors=5")
include_directories(/usr/include/eigen3)
set(main_SOURCES main.cpp)

add_executable(main ${main_SOURCES})
