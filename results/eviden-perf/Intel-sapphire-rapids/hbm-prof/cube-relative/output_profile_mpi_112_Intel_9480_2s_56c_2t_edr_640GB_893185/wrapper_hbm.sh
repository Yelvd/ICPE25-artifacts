#!/bin/bash

rank=$OMPI_COMM_WORLD_LOCAL_RANK
diviseur=$(( rank / 14 + 8 ))
numactl -C $rank --membind $diviseur $@ 

