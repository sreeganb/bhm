#!/usr/bin/bash

rm -rf build/
mkdir build/
cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -DIMP_DIR=/home/sree/work/imp_release/ 
make -j8
make install
