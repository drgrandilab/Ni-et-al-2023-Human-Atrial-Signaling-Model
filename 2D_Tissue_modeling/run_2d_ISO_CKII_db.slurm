#!/bin/bash
#### By Haibo
# A job submission script for running a hybrid MPI/OpenMP job on  
### Job name
#SBATCH --job-name=MPI_TWOD

### Number of compute nodes
#SBATCH --nodes=8

### Number of MPI ranks
##SBATCH --ntasks=1

### Number of MPI ranks per node
#SBATCH --ntasks-per-node=1

### Number of tasks per MPI rank
#SBATCH --cpus-per-task=8

### For how long
#SBATCH --time=3-0:5:0 #  2-0 #unlimited

### Output and error f
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err

###  Assign nodes
##SBATCH --nodelist=c0-[3,8-11]

export OMP_NUM_THREADS=8 #$SLURM_CPUS_PER_TASK
date "+%H:%M:%S   %d/%m/%y"
#srun  --mpi=pmi2 ./Atria_2D_MPI_Ghost Total_time 6000 BCL 1000

# echo ${1} 
# echo ${2} 
# echo ${3}
# ./ONE_D_MPI_Ghost ICs Specific  ISO_con 0.1 Total_time 5600 S1_number 6 S2 330.51
# mv OneD_output.dat.0 OneD_output.dat.0.ISO_con.0.1
# mv cv_out_file.dat cv_out_file.dat.ISO_con.0.1
# mv OneD_output_diff_coef.dat.0 OneD_output_diff_coef.dat.0.ISO_con.0.1
# python2.7 run_VW.py ${1} ${2} ${3}
# python2.7 run_VW_2Hz.py ${1} ${2} ${3}
ISO_con=$1  #0.1
CaMKII=$2  #CaMKII_db
BCL=$3  #1000
# CaMKII_K1_NaV_GapG=No_CaMKII_GapG
folder=BCL.$BCL.ISO.$ISO_con.CaMKII.$CaMKII  #.$CaMKII_K1_NaV_GapG


echo $folder
# mkdir $folder
# ./Atria_2D_MPI_Ghost Total_time 25000 BCL $BCL S2 $BCL S1_number 4 ICs Cell_Specific ISO_con $ISO_con CaMKII $CaMKII Output_Folder $folder 
# # CaMKII_K1_NaV_GapG $CaMKII_K1_NaV_GapG

# mv *.bin $folder


date "+%H:%M:%S   %d/%m/%y"


