#!/bin/bash
#SBATCH -p Intel_9480_2s_56c_2t_edr_640GB
#SBATCH -N 1
#SBATCH -t 6:0:0
#SBATCH -J HEMOCELL
#SBATCH --exclusive
#SBATCH --account=pro_2020_compbiomed2


envfile=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/Intel-sapphire-rapids/env.sh
datadir=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/General
binexe=/home_nfs/bhamitono/projects/compbiomed2/wp4/hemocell/HemoCell-dev-develop/build_intel_sapphirerapids/cube_gnu_sapphirerapids


export JMP=2
export SLOTS=224

#########
# Job info
ntasks=112
omp=1
ppn=112

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


WORKDIR=${PWD}/config-${hemoconfig}/output_energy_mpi_${ntasks}_${SLURM_JOB_PARTITION}_${SLURM_JOBID}

# environment : compiler, mpi, ..
source $envfile

# create workdir
mkdir -p    ${WORKDIR}
cp $envfile ${WORKDIR}/
cp $0       ${WORKDIR}/

# copy data
cp -r ${datadir}/* ${WORKDIR}/
cp config-${hemoconfig}.xml ${WORKDIR}/
cp /scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/Intel-sapphire-rapids/wrapper_hbm.sh ${WORKDIR}/
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
echo "OFF" && clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/turbo_off.sh
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


/home_nfs/bhamitono/software/ipmi-energy-record/ipmi-energy-record.sh -p 0.5 -o ./report_${SLURM_JOBID}.csv -- mpirun -n ${ntasks} --mca pml ucx -x UCX_NET_DEVICES=all --hostfile ./hostfile_$SLURM_JOBID --map-by ppr:${ppn}:node:PE=$JMP --bind-to hwthread ./wrapper_hbm.sh ${binexe} config-${hemoconfig}.xml | tee -a bull_log.out


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
