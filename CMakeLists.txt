cmake_minimum_required(VERSION 3.7)
project(libint-eigen)

set(CMAKE_CXX_STANDARD 11)


# ENABLE LIBINT2
include_directories(/usr/local/libint/2.2.0/include)
include_directories(/usr/local/libint/2.2.0/include/libint2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/opt/local/include/eigen3")


# ENABLE THIS PROJECT
add_executable(libint-eigen
        src/main.cpp
        src/libint-wrapper.cpp
        src/libint-wrapper.hpp)


# ENABLE LIBINT2
target_link_libraries(libint-eigen /usr/local/libint/2.2.0/lib/libint2.a)