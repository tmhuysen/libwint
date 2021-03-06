# libwint v3.0.0

[![Build Status](https://travis-ci.org/GQCG/libwint.svg?branch=master)](https://travis-ci.org/GQCG/libwint)

A C++ library that is a wrapper around libint2. The overlap, kinetic, nuclear and Coulomb repulsion integrals are stored in Eigen3 matrices. For the Coulomb repulsion integrals, the corresponding tensor should be accessed using **chemist's notation**. The library also provides support for integral transformations.


## Dependencies
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-blue.svg)](www.boost.org)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-blue.svg)](https://github.com/evaleev/libint)
[![cpputil Dependency](https://img.shields.io/badge/cpputil-1.2.1+-blue.svg)](https://github.com/GQCG/cpputil)



## Installation
To install this library:
1. clone the master branch, which contains the latest release

        git clone https://github.com/GQCG/libwint.git --branch master --single-branch
        cd libwint

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libwint.a` will be installed in `prefix/libwint/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/libwint/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(libwint x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `libwint_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `libwint` to be used in your `target_link_libraries`.
