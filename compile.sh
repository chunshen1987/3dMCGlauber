#!/usr/bin/env bash

num_of_cores=$1
if [ -z "$num_of_cores" ]
then
    num_of_cores=4
fi

mkdir -p build
cd build
rm -fr *
cmake ..
make -j$num_of_cores
make install
cd ..
rm -fr build
