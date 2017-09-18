# Main CMakeLists (top-level directory)

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(libint-eigen LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)  # LibInt2 needs C++ 11


# Default variables
set(LIBINT_PREFIX /usr/local/libint2/)
set(EIGEN_PREFIX /opt/local/include/)


# source files for the library
add_subdirectory(src)

# source files for the examples
add_subdirectory(examples)


# Set up testing
enable_testing()
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_subdirectory(tests)