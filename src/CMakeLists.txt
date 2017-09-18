# src-level CMakeLists.txt

# We want to make a (static) library called "libint-eigen", based on the given source files
add_library(libint-eigen
        integrals.cpp
        Molecule.cpp
        Basis.cpp)


# Include this project's headers
target_include_directories(libint-eigen PRIVATE ../include)

# Include the libint2 headers
target_include_directories(libint-eigen PRIVATE ${LIBINT_PREFIX}/include)
target_include_directories(libint-eigen PRIVATE ${LIBINT_PREFIX}/include/libint2)

# LibInt2 specifically asks to include eigen3 headers
target_include_directories(libint-eigen PRIVATE ${EIGEN_PREFIX}/eigen3)


# Link against the libint2 library
target_link_libraries(libint-eigen PRIVATE ${LIBINT_PREFIX}/lib/libint2.a)
