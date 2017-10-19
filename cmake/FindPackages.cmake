# In this CMake file, we will find all required packages

# Find the Eigen3 package - needed for this project, and for libint2
# Options to find this package from (http://eigen.tuxfamily.org/dox/TopicCMakeGuide.html)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# Find the boost package - needed for unittests
find_package(Boost REQUIRED)


# There is no find_package(libint2) available, so we'll work with an intermediary solution.
set(libint2_INCLUDE_DIRS ${LIBINT_PREFIX}/include ${LIBINT_PREFIX}/include/libint2)
find_library(libint2 libint2.a HINTS ${LIBINT_PREFIX}/lib)
