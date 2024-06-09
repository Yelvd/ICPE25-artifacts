#!/bin/bash


for conf in "fixed" "relative" "fixed-large"
do
    for j in 16 32 64 80
    do
        for k in {1..1}
        do
            sbatch job_time.sh $conf $j
        done
    done
done

for conf in "fixed" "relative" "fixed-large"
do
    for j in 16 32 64 80
    do
        for k in {1..1}
        do
            sbatch job_energy.sh $conf $j
        done
    done
done

