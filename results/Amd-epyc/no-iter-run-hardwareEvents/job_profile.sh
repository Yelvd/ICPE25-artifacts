#!/bin/bash
#SBATCH -p AMD_9654_2s_96c_2t_edr_768GB_4800
#SBATCH -N 1
#SBATCH -t 6:0:0
#SBATCH -J HEMOCELL
#SBATCH --exclusive
#SBATCH --account=pro_2020_compbiomed2


envfile=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/Amd-epyc/env.sh
datadir=/scratch_lustre_lscoe/users/bhamitono/hemocell/benchmark-setup/General
binexe=/home_nfs/bhamitono/projects/compbiomed2/wp4/hemocell/HemoCell-dev-develop/build_amd_znver4/cube_amd_znver4


export JMP=2
export SLOTS=384

#########
# Job info
ntasks=192
omp=1
ppn=192

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


WORKDIR=${PWD}/config-${hemoconfig}/output_profile_mpi_${ntasks}_${SLURM_JOB_PARTITION}_${SLURM_JOBID}

# environment : compiler, mpi, ..
source $envfile
module use /software/modulefiles/debug-profile/likwid
module load x86_64/5.3.0

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
echo "Perf Event Paranoid to -1"
clush -bw $SLURM_JOB_NODELIST sudo /home_nfs/script/admin/perf_event_paranoid.sh -1
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

mpirun -n ${ntasks} --mca pml ucx -x UCX_NET_DEVICES=all --hostfile ./hostfile_$SLURM_JOBID --map-by ppr:${ppn}:node:PE=$JMP --bind-to hwthread likwid-perfctr -f -g FLOPS_DP -g CACHE -g CLOCK -g DATA -g L2 -g L3 -g MEM1 -g MEM2 -o likwid_prof_rank_%r.txt ${binexe} config-${hemoconfig}.xml | tee -a bull_log.out

#perf stat -e branch-instructions,branch-loads,branch-misses,branch-load-misses,cache-misses,cache-references,cpu-cycles,instructions,stalled-cycles-backend,stalled-cycles-frontend,cpu-clock,cpu-migrations,page-faults,task-clock,L1-dcache-loads,L1-dcache-prefetches,L1-dcache-load-misses,L1-icache-loads,L1-icache-load-misses,fp_arith_inst_retired.128b_packed_double,fp_arith_inst_retired.128b_packed_single,fp_arith_inst_retired.256b_packed_double,fp_arith_inst_retired.256b_packed_single,fp_arith_inst_retired.double,fp_arith_inst_retired.packed,fp_arith_inst_retired.scalar,fp_arith_inst_retired.scalar_double,fp_arith_inst_retired.scalar_single,fp_arith_inst_retired.single,inst_retired.x87,duration_time --cpu 0 -o perf_report.log mpirun -n ${ntasks} --mca pml ucx -x UCX_NET_DEVICES=all --hostfile ./hostfile_$SLURM_JOBID --map-by ppr:${ppn}:node:PE=$JMP --bind-to hwthread ${binexe} config-${hemoconfig}.xml | tee -a bull_log.out


#perf stat -e branch-instructions,branch-loads,branch-misses,branch-load-misses,cache-misses,cache-references,cpu-cycles,instructions,stalled-cycles-backend,stalled-cycles-frontend,cpu-clock,cpu-migrations,page-faults,task-clock,L1-dcache-loads,L1-dcache-prefetches,L1-dcache-load-misses,L1-icache-loads,L1-icache-load-misses,duration_time,fp_ops_retired_by_type.scalar_all,fp_ops_retired_by_type.vector_all,fp_ops_retired_by_width.pack_128_uops_retired,fp_ops_retired_by_width.pack_256_uops_retired,fp_ops_retired_by_width.pack_512_uops_retired,fp_ops_retired_by_width.scalar_uops_retired,fp_ops_retired_by_width.x87_uops_retired,fp_pack_ops_retired.fp128_all,fp_pack_ops_retired.fp256_all,fp_ret_x87_fp_ops.all,duration_time --cpu 0 -o perf_report.log mpirun -n ${ntasks} --mca pml ucx -x UCX_NET_DEVICES=all --hostfile ./hostfile_$SLURM_JOBID --map-by ppr:${ppn}:node:PE=$JMP --bind-to hwthread ${binexe} config-${hemoconfig}.xml | tee -a bull_log.out

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
