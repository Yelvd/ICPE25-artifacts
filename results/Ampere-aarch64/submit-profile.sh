#!/bin/bash


for conf in "fixed" "relative" "fixed-large"
do
    for j in 64 80
    do
        for k in {1..2}
        do
            sbatch job_profile.sh $conf $j
        done
    done
done


