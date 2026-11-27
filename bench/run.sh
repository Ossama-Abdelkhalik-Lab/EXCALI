#!/bin/bash
# Build every variant with the given compiler and flags, then time them in
# interleaved rounds on generated cases. Output: CSV lines
# compiler,flags,variant,round,set,single_ns,batch_ns
MODFLAG=${MODFLAG:--J}; FC=${FC:-gfortran}; TAG=${TAG:-base}; FLAGS=${FLAGS:-"-O3 -ffree-line-length-none -fno-fast-math"}
VARIANTS=${VARIANTS:-$(ls bench/variants)}
ROUNDS=${ROUNDS:-3}
out=bench/out/$TAG; mkdir -p $out
for v in $VARIANTS; do
  d=bench/variants/$v; o=$out/obj_$v; mkdir -p $o
  if $FC $FLAGS $MODFLAG $o -I $d -c $d/excali_solver.f90 -o $o/ex.o 2> $o/build.log && \
     $FC $FLAGS -I $o bench/tgen.f90 $o/ex.o -o $out/tg_$v 2>> $o/build.log; then :; else echo "BUILD FAILED $v"; cat $o/build.log | head -5; fi
done
for r in $(seq 1 $ROUNDS); do
  for v in $VARIANTS; do
    [ -x $out/tg_$v ] || continue
    $out/tg_$v | awk -v fc="$FC" -v tag="$TAG" -v v=$v -v r=$r '{printf "%s,%s,%s,%s,%s,%s,%s\n",fc,tag,v,r,$1,$3,$5}'
  done
done
