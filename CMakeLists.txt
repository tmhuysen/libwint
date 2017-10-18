# top-level CMakeLists

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(libint-wrapper LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)


# Default variables decision trees
#   LIBINT_PREFIX
if(NOT LIBINT_PREFIX)
    if(DEFINED ENV{LIBINT_PREFIX})
        set(LIBINT_PREFIX $ENV{LIBINT_PREFIX})
    else()
        message(WARNING "You did not specify -DLIBINT_PREFIX, nor is LIBINT_PREFIX set in your environment. Proceeding with default value for LIBINT_PREFIX.")
        set(LIBINT_PREFIX /usr/local/libint/2.3.1)
    endif()
endif()

#   INSTALLATION_PREFIX
if(NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# Find the Eigen3 package (http://eigen.tuxfamily.org/dox/TopicCMakeGuide.html)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# source files for the library
add_subdirectory(src)


# Set up testing
enable_testing()
find_package(Boost REQUIRED)
add_subdirectory(tests)
