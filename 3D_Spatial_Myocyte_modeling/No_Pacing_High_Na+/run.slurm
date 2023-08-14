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
ISO=$3
mkdir -p global_result # create folders for data ouput
mkdir -p BinaryFiles
time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/AF.dat.BCL$BCL.AF.1.ISO.0.CKII_ih.0.CKII_db.1" $ISO
mkdir -p BCL.$BCL.CKII_db
rsync -a steady_state_output steady_state_init global_result BinaryFiles AF.BCL.$BCL.CKII_db.ISO.$ISO



# mkdir -p global_result # create folders for data ouput
# mkdir -p BinaryFiles
# time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/Data_BCL.$BCL.ISO.0.CaMKII_inh.0.CaMKII_db.0_para.dat"

# mkdir -p BCL.$BCL.Normal
# rsync -a steady_state_output steady_state_init global_result BinaryFiles BCL.$BCL.Normal



mkdir -p global_result # create folders for data ouput
mkdir -p BinaryFiles
time ./pace_2 $BCL $Cycles "pool_tubule/tub_input_ver2_20.txt" "CaMKII_levels/CaMKII_inhbition.dat" $ISO
mkdir -p BCL.$BCL.CKII_inh
rsync -a steady_state_output steady_state_init global_result BinaryFiles AF.CL.$BCL.CKII_inh.ISO.$ISO
}


run 250 112 0
run 250 112 1
# run 333.333 84
# run 500 56
# run 1000 28

