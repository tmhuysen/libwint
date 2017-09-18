# libint-eigen

A C++ library that stores LibInt2 calculated integrals in eigen3 matrices

## dependencies
[![libint Dependency](https://img.shields.io/badge/libint-2.2.0+-blue.svg)](https://github.com/evaleev/libint)
[![eigen3 Dependency](https://img.shields.io/badge/eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## installation

To install this library, clone this repository

    git clone git@github.com:lelemmen/libint-eigen.git

and perform an out-of-source build:

    mkdir build && cd build
    cmake -DLIBINT_PREFIX=your_libint_prefix -DEIGEN_PREFIX=your_eigen_prefix ..

where your_libint_prefix is the prefix to your libint2 installation folder, and your_eigen_prefix is the prefix to your eigen3 headers.