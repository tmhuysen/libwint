# examples-level CMakeLists

# Create the example executable
add_executable(example_integrals example_integrals.cpp)


# Obviously, it needs to include this project's headers
target_include_directories(example_integrals PRIVATE ../include)

# Include the eigen3 headers
target_include_directories(example_integrals PRIVATE /opt/local/include/eigen3)

# Include the libint2 headers
target_include_directories(example_integrals PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(example_integrals PRIVATE /usr/local/libint/2.2.0/include/libint2)


# Link with our own library
target_link_libraries(example_integrals libint-eigen)
