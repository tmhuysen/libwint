# Main CMakeLists (top-level directory)

cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(libint-eigen LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)  # LibInt2 needs C++ 11


# SCOPE: We're developing a library that interfaces LibInt2 with eigen3, which has its source file(s) in src
#   For development, however, I like adding a target based on main.cpp, so we will have to compile the source files and then turn them into a library
#   Luckily, cmake does this for us
add_subdirectory(src)  # Note that this subdirectory has its own CMakeLists.txt specifying the target 'libint-eigen-wrapper', which will be our library


# Create the development executable, based on main.cpp
add_executable(libint-eigen main.cpp)

# Obviously, it needs to include this project's headers
target_include_directories(libint-eigen PRIVATE include)

# We will have to include the headers to LibInt2 and eigen3
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include/libint2)
target_include_directories(libint-eigen PRIVATE /opt/local/include/eigen3)

# Furthermore, we have to link with our own library (which is specified in src)
target_link_libraries(libint-eigen libint-eigen-wrapper)

# Q: Why don't we link with the LibInt2 library here?
# A: We have already linked our own library with it, and main.cpp doesn't specify any LibInt2 computations


enable_testing()                # enables testing for this directory and below

add_test(NAME first_test COMMAND first_test_command)