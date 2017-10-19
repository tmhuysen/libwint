# src-level CMakeLists.txt


# We want to make a static library target called "libwrp", based on the project source files
add_library(${LIBRARY_NAME} ${LIBRARY_TYPE} ${PROJECT_SOURCE_FILES})
set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")     # the final library name should be libwrp.a, and not liblibwrp.a


# Configure the library: include headers and link to libraries
include(${CMAKE_SOURCE_DIR}/cmake/ConfigureLibrary.cmake)


# Install the library
include(${CMAKE_SOURCE_DIR}/cmake/InstallLibrary.cmake)
