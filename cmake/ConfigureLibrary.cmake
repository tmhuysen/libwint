# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(libint-wrapper PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the libint2 headers
target_include_directories(libint-wrapper PRIVATE ${LIBINT_PREFIX}/include)
target_include_directories(libint-wrapper PRIVATE ${LIBINT_PREFIX}/include/libint2)

# Link against the libint2 library
target_link_libraries(libint-wrapper PUBLIC ${LIBINT_PREFIX}/lib/libint2.a)

# Include the Eigen3 headers
target_link_libraries(libint-wrapper PUBLIC Eigen3::Eigen)