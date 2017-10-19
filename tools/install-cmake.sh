#!/usr/bin/env bash

# Travis CI's CMake version is only 3.2.2.
# This is a workaround that upgrades to a newer version of CMake.
# Credits to (https://github.com/travis-ci/travis-ci/issues/7437)


cd
mkdir ${HOME}/usr
export PATH="$HOME/usr/bin:$PATH"
wget https://cmake.org/files/v3.9/cmake-3.9.1-Linux-x86_64.sh
chmod +x cmake-3.9.1-Linux-x86_64.sh
./cmake-3.9.1-Linux-x86_64.sh --prefix=$HOME/usr --exclude-subdir --skip-license
cd
