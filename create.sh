#!/bin/bash

rm -rf build/
mkdir build/
cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -DIMP_DIR=/opt/homebrew/Cellar/imp/2.20.2_1/lib/cmake/IMP/ 
make -j8
make install
