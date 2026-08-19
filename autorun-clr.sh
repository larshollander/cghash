#!/bin/bash

REPS=5
N_INPUTS="1000"
N_DELTAS="10"
D_LENGTHS=("2048" "3800" "6800")
A_LENGTHS=("266" "274" "280" "286" "292" "297" "302" "307" "311")

cwd=$(pwd)

if [ -e "results-clr.txt" ]; then
  rm -i results-clr.txt
fi

for _ in $(seq 1 $REPS);
  do for d_length in "${D_LENGTHS[@]}";
    do for i in {1..8};
      do $cwd/c/a.out $N_INPUTS $N_DELTAS $d_length ${A_LENGTHS[$i]} crl$i csld$i csmld$i >> results-clr.txt
    done
  done
done

python plot_runtimes.py results-clr.txt
