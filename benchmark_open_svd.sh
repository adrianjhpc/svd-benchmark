#!/bin/bash --login
#
#PBS -N svd_open_bench
#PBS -l select=1:ncpus=36
# Parallel jobs should always specify exclusive node access
#PBS -l place=excl
#PBS -l walltime=1:00:0
#PBS -A ec110-guest

module load intel-compilers-17
module load intel-mpi-17
module load intel-cmkl-17

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/sc004/adrianjc/scalapack/2.1.0/lib

export OMP_NUM_THREADS=8
cd /lustre/home/sc004/adrianjc/svd

for i in 16 32 64 128 256 512 1024 2048 4096
do
  mpirun -n 4 ./scalapack_open $i $i 2 2
  mpirun -n 4 ./scalapack_open $i $i 2 2
  mpirun -n 4 ./scalapack_open $i $i 2 2
  mpirun -n 4 ./scalapack_open $i $i 4 1
  mpirun -n 4 ./scalapack_open $i $i 4 1
  mpirun -n 4 ./scalapack_open $i $i 4 1
  mpirun -n 8 ./scalapack_open $i $i 4 2
  mpirun -n 8 ./scalapack_open $i $i 4 2
  mpirun -n 8 ./scalapack_open $i $i 4 2
  mpirun -n 8 ./scalapack_open $i $i 8 1
  mpirun -n 8 ./scalapack_open $i $i 8 1
  mpirun -n 8 ./scalapack_open $i $i 8 1
done


#for i in 1024 2048 4096 8192 16384 32768
#do
 # mpirun -n 8 ./scalapack_open $i $(( $i/8 ))
 # mpirun -n 8 ./scalapack_open $i $(( $i/8 ))
 # mpirun -n 8 ./scalapack_open $i $(( $i/8 ))
#done

