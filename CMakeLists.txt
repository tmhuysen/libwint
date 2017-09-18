# top-level CMakeLists

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(libint-eigen LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)


# Default variables decision trees
#   libint
if(NOT LIBINT_PREFIX)
    if(DEFINED ENV{LIBINT_PREFIX})
        set(LIBINT_PREFIX $ENV{LIBINT_PREFIX})
    else()
        message(WARNING "You did not specify -DLIBINT_PREFIX, nor is LIBINT_PREFIX set in your environment. Proceeding with default value for LIBINT_PREFIX.")
        set(LIBINT_PREFIX /usr/local/libint/2.2.0)
    endif()
endif()

#   eigen
if(NOT EIGEN_PREFIX)
    if(DEFINED ENV{EIGEN_PREFIX})
        set(EIGEN_PREFIX $ENV{EIGEN_PREFIX})
    else()
        message(WARNING "You did not specify -DEIGEN_PREFIX, nor is EIGEN_PREFIX set in your environment. Proceeding with default value for EIGEN_PREFIX.")
        set(EIGEN_PREFIX /opt/local/include)
    endif()
endif()


# source files for the library
add_subdirectory(src)

# source files for the examples
add_subdirectory(examples)


# Set up testing
enable_testing()
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_subdirectory(tests)



#