#!/bin/bash


for conf in "fixed" "relative" "fixed-large"
do
    for j in 56 112
    do
        for k in {1..1}
        do
            sbatch job_profile.sh $conf $j
        done
    done
done


