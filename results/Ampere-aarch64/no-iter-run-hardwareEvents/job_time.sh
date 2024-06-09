#!/bin/bash
#SBATCH -p Ampere_Q8030_1s_80c_1t_edr_256GB_3200
#SBATCH -N 1
#SBATCH -t 6:0:0
#SBATCH -J HEMOCELL
#SBATCH --exclusive
#SBATCH --account=pro_2020_compbiomed2
#SBATCH --reservation=sbatch_test


envfile=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/Ampere-aarch64/env.sh
datadir=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/General
binexe=/home_nfs/bhamitono/projects/compbiomed2/wp4/hemocell/HemoCell-dev-develop/build_ampere_gnu/cube_aarch64


export JMP=1
export SLOTS=80

#########
# Job info
ntasks=80
omp=1
ppn=80

if [ $# -ge 1 ]
then
    hemoconfig=$1
    ntasks=$2
fi

if [ "$hemoconfig" = "fixed" ]; then
    echo "CONFIG = FIXED"
elif [ "$hemoconfig" = "relative" ];then
    echo "CONFIG = RELATIVE"
elif [ "$hemoconfig" = "fixed-large" ]; then
    echo "CONFIG = FIXED-LARGE"
else
    echo "=========================================="
    echo "Something went wrong ."
    echo "=========================================="
    exit 1
fi


WORKDIR=${PWD}/config-${hemoconfig}/output_time_mpi_${ntasks}_${SLURM_JOB_PARTITION}_${SLURM_JOBID}

# environment : compiler, mpi, ..
source $envfile

# create workdir
mkdir -p    ${WORKDIR}
cp $envfile ${WORKDIR}/
cp $0       ${WORKDIR}/

# copy data
cp -r ${datadir}/* ${WORKDIR}/
cp config-${hemoconfig}.xml ${WORKDIR}/
# move to workdir
cd ${WORKDIR}

#########
# CPU settings: to be commented if outside atos
echo "Zone reclaim"
clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/zone_reclaim 1
echo "Numa balancing"
clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/enable_numa_balancing
echo "Transparent HP"
clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/enable_transparent_hugepages_rhel7
echo "Drop Cache"
clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/dropcache
echo "Turbo"
#echo "ON" && clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/turbo_on.sh
#echo "OFF" && clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/turbo_off.sh
#FREQ=$1
#echo "Frequency set" $FREQ
#clush -bw $SLURM_JOB_NODELIST "sudo /home_nfs/script/admin/pwrb -f $FREQ"


ulimit -s unlimited
#########

################
# OPENMP
export OMP_NUM_THREADS=${omp}
if [ "$OMP_NUM_THREADS" -ge 2 ]
then
    export tmp=$(( OMP_NUM_THREADS -1 ))
    export GOMP_AFFINITY=$(seq -s, 0 1 $tmp )
fi
################


################

################
# Launch
nodeset -e $SLURM_NODELIST | tr ' ' '\n' > ./hostfilee
# slots refers to the number of hyperthreads and not number of tasks
while read -r line ; do echo "$line slots=$SLOTS"; done < hostfilee > hostfile_$SLURM_JOBID
rm hostfilee


echo ' #####################################'
echo ' Begin : '
echo ' #####################################'

/usr/bin/time -f %e -o bull_time.out mpirun -n ${ntasks} --mca pml ucx -x UCX_NET_DEVICES=all --hostfile ./hostfile_$SLURM_JOBID --map-by ppr:${ppn}:node:PE=$JMP --bind-to hwthread ${binexe} config-${hemoconfig}.xml | tee -a bull_log.out


# clean
rm RBC.xml
rm RBC.pos
rm PLT.xml
rm PLT.pos
rm -rf tmp_1 

vartime=`cat bull_time.out`
################

echo "BULL MPI TASKS   = ${ntasks}"
echo "BULL Bin name    = ${binexe} "
echo "BULL config      = $hemoconfig "
echo "BULL PARTITION   = ${SLURM_JOB_PARTITION}"
echo "BULL Time        = ${vartime}"

cp ${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOBID}.out .

################
