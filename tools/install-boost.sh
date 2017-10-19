#!/usr/bin/env bash

# We need relatively recent boost libraries, and APT doesn't provide these.
# This script installs the header files in /usr/local and exports BOOST_ROOT.


cd
mkdir /tmp/boost && cd /tmp/boost
curl -OL "https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2"
cd /usr/local/
sudo tar --bzip2 -xf /tmp/boost/boost_1_65_1.tar.bz2
export BOOST_ROOT=/usr/local/boost_1_65_1
cd
