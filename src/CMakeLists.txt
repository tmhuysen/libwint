# src-level CMakeLists.txt

# We want to make a (static) library called "libint-eigen", based on the given source files
add_library(libint-eigen
        integrals.cpp
        Molecule.cpp
        Basis.cpp)


# Include this project's headers
target_include_directories(libint-eigen PRIVATE ../include)

# Include the libint2 headers
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(libint-eigen PRIVATE /usr/local/libint/2.2.0/include/libint2)

# LibInt2 specifically asks to include eigen3 headers
target_include_directories(libint-eigen PRIVATE /opt/local/include/eigen3)


# Link against the libint2 library
target_link_libraries(libint-eigen /usr/local/libint/2.2.0/lib/libint2.a)
