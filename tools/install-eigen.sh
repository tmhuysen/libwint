#!/usr/bin/env bash

# We need a recent version of Eigen3, but this is not available through APT.
# This script installs the latest version of Eigen3.

cd
hg clone https://bitbucket.org/eigen/eigen /tmp/eigen
mkdir /tmp/build-eigen && cd /tmp/build-eigen
cmake . /tmp/eigen
make && sudo make install
cd
