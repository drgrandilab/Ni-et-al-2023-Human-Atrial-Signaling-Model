#!/bin/bash
#SBATCH --job-name=T10_1hz
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
#SBATCH --cpus-per-task=56
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

make clean; make pace_2; 

run() {
BCL=$1
Cycles=$2
mkdir -p global_result # create folders for data ouput
mkdir -p BinaryFiles
time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/test_1.dat"

mkdir -p BCL.$BCL.Normal
rsync -a steady_state_output steady_state_init global_result BinaryFiles BCL.$BCL.test_1

mkdir -p global_result # create folders for data ouput
mkdir -p BinaryFiles
time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/test_2.dat"
mkdir -p BCL.$BCL.CKII_db
rsync -a steady_state_output steady_state_init global_result BinaryFiles BCL.$BCL.test_2


# mkdir -p global_result # create folders for data ouput
# mkdir -p BinaryFiles
# time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/CaMKII_inhbition.dat"
# mkdir -p BCL.$BCL.CKII_inh
# rsync -a steady_state_output steady_state_init global_result BinaryFiles BCL.$BCL.CKII_inh
}



run 250 112
# run 333.333 84
# run 500 56
# run 1000 28

