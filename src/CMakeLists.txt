# src-level CMakeLists.txt


# We want to make a static library target called "libint-wrapper", based on the project source files
add_library(${LIBRARY_NAME} ${LIBRARY_TYPE} ${PROJECT_SOURCE_FILES})


# Configure the library: include headers and link to libraries
include(${CMAKE_SOURCE_DIR}/cmake/ConfigureLibrary.cmake)


# Install the library
include(${CMAKE_SOURCE_DIR}/cmake/InstallLibrary.cmake)

message(WARNING ${PROJECT_NAME_UPPERCASE})
message(WARNING ${PROJECT_VERSION_MAJOR})