# src-level CMakeLists.txt

# We want to make a (static) library called "libint-wrapper", based on the given source files
add_library(libint-wrapper
        integrals.cpp
        Molecule.cpp
        Basis.cpp)


# Include this project's headers
target_include_directories(libint-wrapper PRIVATE ../include)

# Include the libint2 headers
target_include_directories(libint-wrapper PRIVATE ${LIBINT_PREFIX}/include)
target_include_directories(libint-wrapper PRIVATE ${LIBINT_PREFIX}/include/libint2)

# Link against the libint2 library
target_link_libraries(libint-wrapper PUBLIC ${LIBINT_PREFIX}/lib/libint2.a)

# Include the Eigen3 headers
target_link_libraries(libint-wrapper PUBLIC Eigen3::Eigen)


# Installation of the target library and its associated header files
install(TARGETS libint-wrapper DESTINATION ${INSTALLATION_PREFIX}/lib)
install(DIRECTORY ../include/ DESTINATION ${INSTALLATION_PREFIX}/include/libint-wrapper)