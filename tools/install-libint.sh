#!/usr/bin/env bash

# We can't install libint2 (2.2.0+) through APT.
# This script installs a recent version of libint2, carefully following the steps from the libint2 wiki (https://github.com/evaleev/libint/wiki)

# Installing libint consists of two parts: compiling the libint compiler and compiling the libint library.
#   1. Compile the libint compiler
cd
mkdir /tmp/libint && cd /tmp/libint
curl -OL "https://github.com/evaleev/libint/archive/v2.3.1.tar.gz"
tar -xvzf v2.3.1.tar.gz
cd libint-2.3.1
./autogen.sh
mkdir build && cd build
../configure CXXFLAGS=-I${BOOST_ROOT}   # the environment variable BOOST_ROOT was set in the install-boost.sh script

#   2. Generate and compile the libint library
make export
tar -xvzf libint-2.3.1.tgz
cd libint-2.3.1
./configure CXXFLAGS=-I${BOOST_ROOT}    # the environment variable BOOST_ROOT was set in the install-boost.sh script
make
sudo make install

# Finally, since libint often complains about reading the basis sets, specify a LIBINT_DATA_PATH to help libint find the files.
export LIBINT_DATA_PATH=/usr/local/libint/2.3.1/share/libint/2.3.1/basis
cd
