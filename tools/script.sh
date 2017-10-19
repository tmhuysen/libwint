#!/usr/bin/env bash

# This is the build script that determines if the build fails or passes.

cd
cd /home/travis/build/lelemmen/libint-wrapper
mkdir build && cd build
cmake .. -DBOOST_ROOT=${BOOST_ROOT}   # the environment variable BOOST_ROOT was set in the install-boost.sh script
make && make test && sudo make install
cd
