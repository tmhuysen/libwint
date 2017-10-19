# libwrp

[![Build Status](https://travis-ci.org/lelemmen/libwrp.svg?branch=master)](https://travis-ci.org/lelemmen/libwrp)

A C++ library that stores libint2 calculated overlap, kinetic, nuclear and Coulomb repulsion integrals in Eigen3 matrices. For the Coulomb repulsion integrals, the corresponding tensor should be accessed using **chemist's notation**.


## Dependencies
[![libint2 Dependency](https://img.shields.io/badge/LibInt-2.3.1+-blue.svg)](https://github.com/evaleev/libint)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## Installation
To install this library:
1. clone this repository

        git clone git@github.com:lelemmen/libwrp.git


2. perform an out-of-source build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=your_wanted_installation_prefix -DLIBINT_PREFIX=your_libint_prefix ..
        make && make test

    where
    * `your_wanted_installation_prefix` is the installation prefix you want the library (installed at `your_wanted_installation_prefix/lib`) and its headers (`your_wanted_installation_prefix/include`) to be installed. `your_wanted_installation_prefix` defaults to `/usr/local`;
    * `your_libint_prefix` is the prefix to your libint2 installation folder, which defaults to `${LIBINT_PREFIX}` if available in `env`, and if not available, defaults to `/usr/local/libint/2.3.1`;


3. finally, install (possibly with admin privileges via `sudo`)

        make install


## Usage
Basic use of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file

    find_package(libwrp x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `libwrp_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `libwrp` to be used in your `target_link_libraries`.
