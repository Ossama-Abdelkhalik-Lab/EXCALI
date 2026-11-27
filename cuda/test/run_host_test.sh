#!/bin/bash
# Host check of the CUDA port against Fortran excali_one on truth CSVs.
# usage: cuda/test/run_host_test.sh <csv>... (from the repository root)
# Builds both with and without FMA contraction; without it the two must be bit-identical.
set -e
B=${BUILD:-build/cuda_host}; mkdir -p $B/obj $B/objn
FC=${FC:-gfortran}; CXX=${CXX:-g++}
$FC -O3 -ffree-line-length-none -fno-fast-math -J $B/obj -c src/excali_solver.f90 -o $B/obj/ex.o
$FC -O2 -I $B/obj cuda/test/fortran_driver.f90 $B/obj/ex.o -o $B/fortran_driver
$FC -O3 -ffp-contract=off -ffree-line-length-none -fno-fast-math -J $B/objn -c src/excali_solver.f90 -o $B/objn/ex.o
$FC -O2 -ffp-contract=off -I $B/objn cuda/test/fortran_driver.f90 $B/objn/ex.o -o $B/fortran_driver_nofma
$CXX -O3 -std=c++17 cuda/test/host_driver.cpp -o $B/host_driver
$CXX -O3 -std=c++17 -ffp-contract=off cuda/test/host_driver.cpp -o $B/host_driver_nofma
for f in "$@"; do
  n=$(basename $f .csv)
  $B/fortran_driver $f $B/f_$n.txt; $B/host_driver $f $B/c_$n.txt
  $B/fortran_driver_nofma $f $B/fn_$n.txt; $B/host_driver_nofma $f $B/cn_$n.txt
  echo "== $n, default contraction"; python3 cuda/test/compare.py $B/c_$n.txt $B/f_$n.txt
  echo "== $n, -ffp-contract=off (must be identical)"; python3 cuda/test/compare.py $B/cn_$n.txt $B/fn_$n.txt
done
