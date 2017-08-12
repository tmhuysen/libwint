cmake_minimum_required(VERSION 3.7)
project(libint_tutorial)

set(CMAKE_CXX_STANDARD 11)  # C++11 standard library


# ENABLE LIBINT2
include_directories(/usr/local/libint/2.2.0/include)
include_directories(/usr/local/libint/2.2.0/include/libint2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/opt/local/include/eigen3")
#   Make sure to set Run > Edit Configurations > Environment variables
#       LIBINT_DATA_PATH    /usr/local/libint/2.2.0/include/libint2/basis


set(SOURCE_FILES main.cpp)
add_executable(libint_tutorial ${SOURCE_FILES})  # Makes an executable, based on the given source files


# ENABLE LIBINT2
target_link_libraries(libint_tutorial /usr/local/libint/2.2.0/lib/libint2.a)