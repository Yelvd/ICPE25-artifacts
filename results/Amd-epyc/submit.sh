#!/bin/bash


for conf in "fixed" "relative" "fixed-large"
do
    for j in 32 64 128 256
    do
        for k in {1..4}
        do
            sbatch job_time.sh $conf $j
        done
    done
done

for conf in "fixed" "relative" "fixed-large"
do
    for j in 32 64 128 256
    do
        for k in {1..4}
        do
            sbatch job_energy.sh $conf $j
        done
    done
done

