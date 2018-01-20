# libwrp

[![Build Status](https://travis-ci.org/lelemmen/libwrp.svg?branch=master)](https://travis-ci.org/lelemmen/libwrp)

A C++ library that stores libint2 calculated overlap, kinetic, nuclear and Coulomb repulsion integrals in Eigen3 matrices. For the Coulomb repulsion integrals, the corresponding tensor should be accessed using **chemist's notation**.


## Dependencies
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-blue.svg)](https://github.com/evaleev/libint)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## Installation
To install this library:
1. clone the master branch

        git clone https://github.com/GQCG/libwrp.git
        cd libwrp

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libwrp.a` will be installed in `prefix/libwrp/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/libwrp/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(libwrp x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `libwrp_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `libwrp` to be used in your `target_link_libraries`.
