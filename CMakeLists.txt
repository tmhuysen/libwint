# Main CMakeLists (top-level directory)
cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
project(libint-eigen LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)



enable_testing()                    # enables testing for this directory and below



add_executable(libint-eigen main.cpp)  # create a new target 'libint-eigen'

target_include_directories(libint-eigen PRIVATE include)  # for the target 'libint-eigen',

# add the 'include' directory to those the compiler uses to search for include files
add_subdirectory(src)               # adds the subdirectory 'src' to the build


# add a target with name 'libint-eigen', built from the listed source files
# Include these directories for both main.cpp and libint-wrapper.cpp
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include/libint2)        # FIXME: test if we can leave this open
target_include_directories(libint-eigen PRIVATE /opt/local/include/eigen3)
target_link_libraries(libint-eigen libint-eigen-wrapper)


add_test(NAME first_test COMMAND first_test_command)