#!/bin/bash

for (( s = 1; s <= 4; s++ ))
do \
  for (( n = 1; n <= 30; n++ ))
  do \
    for (( m = 3; m <= n; m +=3))
    do \
      for (( p = 1; p <= 4; p ++ ))
      do \

          mpirun -n $p ./a.out $n $m 0 $s

        done
    done
  done
done
