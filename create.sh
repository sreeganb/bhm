#!/bin/bash

rm -rf build/
mkdir build/
cd build/
# cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -DIMP_DIR=/opt/homebrew/Cellar/imp/2.21.0_5/lib/cmake/IMP/ -DCMAKE_CXX_FLAGS='-std=c++17 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION' -G Ninja 
cmake ../ -DCMAKE_BUILD_TYPE=Release -DIMP_DIR=/opt/homebrew/Cellar/imp/2.21.0_5/lib/cmake/IMP/ -DCMAKE_CXX_FLAGS='-std=c++17 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION' -G Ninja 
make -j8
make install
