#!/bin/bash

REPS=5
N_INPUTS="1000"
N_DELTAS="10"
D_LENGTHS=("2048" "3800" "6800")
A_LENGTHS=("16" "32" "64" "128" "266")
ALGS="sr srd ss ssd ssm ssmd"

cwd=$(pwd)

if [ -e "results-sf.txt" ]; then
  rm -i results-sf.txt
fi

for _ in $(seq 1 $REPS);
  do for d_length in "${D_LENGTHS[@]}";
    do for a_length in "${A_LENGTHS[@]}";
      do $cwd/c/a.out $N_INPUTS $N_DELTAS $d_length $a_length $ALGS >> results-sf.txt
    done
  done
done

python plot_runtimes.py results-sf.txt
