# Main CMakeLists (top-level directory)
cmake_minimum_required(VERSION 3.7)

project(libint-eigen)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)


enable_testing()
add_subdirectory(tests)