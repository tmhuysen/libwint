# examples-level CMakeLists

# Create the example executable
add_executable(example_integrals example_integrals.cpp)


# Obviously, it needs to include this project's headers
target_include_directories(example_integrals PRIVATE ../include)

# Include the eigen3 headers
target_include_directories(example_integrals PRIVATE ${EIGEN_PREFIX}/eigen3)

# Include the libint2 headers
target_include_directories(example_integrals PRIVATE ${LIBINT_PREFIX}/include)
target_include_directories(example_integrals PRIVATE ${LIBINT_PREFIX}/include/libint2)


# Link with our own library
target_link_libraries(example_integrals libint-eigen)
