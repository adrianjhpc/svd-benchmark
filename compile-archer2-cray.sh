#!/bin/bash

CC -fopenmp -O3 -std=c++11 -Wall -I/work/z19/z19/adrianj/compress/eigen/3.3.9/include/eigen3 -o scalapack_libsci scalapack.cc  -lsci_cray_mp -lpthread -lm

echo " "
echo "Compiled libsci version"
echo " "
