# libint-eigen

A C++ library that stores LibInt2 calculated integrals in eigen3 matrices

## dependencies
[![libint Dependency](https://img.shields.io/badge/libint-2.2.0+-blue.svg)](https://github.com/evaleev/libint)
[![eigen3 Dependency](https://img.shields.io/badge/eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## installation

To install this library:
1. clone this repository

        git clone git@github.com:lelemmen/libint-eigen.git


2. perform an out-of-source build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=your_wanted_installation_prefix -DLIBINT_PREFIX=your_libint_prefix -DEIGEN_PREFIX=your_eigen_prefix ..

    where
    * `your_wanted_installation_prefix` is the installation prefix you want the library (installed at `your_wanted_installation_prefix/lib`) and its headers (`your_wanted_installation_prefix/include`) to be installed. `your_wanted_installation_prefix` defaults to `/usr/local`;
    * `your_libint_prefix` is the prefix to your libint2 installation folder, which defaults to `${LIBINT_PREFIX}` if available in `env`, and if not available, defaults to `/usr/local/libint/2.2.0`;
    * `your_eigen_prefix` is the prefix to your eigen3 headers, which defaults to `${EIGEN_PREFIX}` if available in `env`, and if not available, defaults to `/opt/local/include`.


3. finally, install it

        make install