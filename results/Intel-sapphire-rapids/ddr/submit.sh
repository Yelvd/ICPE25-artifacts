#!/bin/bash


for conf in "fixed" "relative" "fixed-large"
do
    for j in 32 56 64 112
    do
        for k in {1..5}
        do
            sbatch job_time.sh $conf $j
        done
    done
done

for conf in "fixed" "relative" "fixed-large"
do
    for j in 32 56 64 112
    do
        for k in {1..5}
        do
            sbatch job_energy.sh $conf $j
        done
    done
done

