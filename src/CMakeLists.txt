# src-level CMakeLists.txt

# We want to make a (static) library called "libint-eigen-wrapper", based on the given source files
add_library(libint-eigen-wrapper
        libint-wrapper.cpp
        Molecule.cpp)


# To make this library, we have to include several headers
#   1. my own header file
target_include_directories(libint-eigen-wrapper PRIVATE ../include)

#   2. LibInt2's header files
target_include_directories(libint-eigen-wrapper PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(libint-eigen-wrapper PRIVATE /usr/local/libint/2.2.0/include/libint2)

#   3. LibInt2 specifically asks to include eigen3 headers, which is convenient as the goal is to store integrals in eigen3 matrices
target_include_directories(libint-eigen-wrapper PRIVATE /opt/local/include/eigen3)


# We also have to link our library with LibInt2's library
target_link_libraries(libint-eigen-wrapper /usr/local/libint/2.2.0/lib/libint2.a)