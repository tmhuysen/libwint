# top-level CMakeLists
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# Specify the project, its version and the used language
project(libint-wrapper
        VERSION 1.0.0
        LANGUAGES CXX)

# To enable the C++11 standard, we have to add flags like -std=gnu+11 or so. CMake can take care of this:
set(CMAKE_CXX_STANDARD 11)



# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

# Find all required packages
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)



# Add the sources for the library
add_subdirectory(src)



# Set up testing
enable_testing()
add_subdirectory(tests)
