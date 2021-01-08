#!/bin/bash
export I_MPI_CXX=icpc

mpicxx -O3 -std=c++11 -Wall -DMKL -I/lustre/home/sc004/adrianjc/puri-psi-dev-build/external/include/eigen3 -I${MKLROOT}/include -o scalapack_mkl scalapack.cc -L${MKLROOT}/lib/intel64 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

echo " "
echo "Compiled MKL version"
echo " "

mpicxx -O3 -std=c++11 -Wall  -I/lustre/home/sc004/adrianjc/puri-psi-dev-build/external/include/eigen3 -o scalapack_open scalapack.cc -L/home/sc004/adrianjc/scalapack/2.1.0/lib/ -lscalapack -L${MKLROOT}/lib/intel64 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm 

echo " "
echo "Compiled open source version"
echo " "
