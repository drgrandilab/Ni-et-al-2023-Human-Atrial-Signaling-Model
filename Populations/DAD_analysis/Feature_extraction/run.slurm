#!/bin/bash
#SBATCH --job-name=EAD_DAD_Analysis
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBTACH --mail-user=haibo.ni02@gmail.com
#SBTACH --mail-type=ALL
#SBTACH --begin=now #+20hour


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK




matlab -nodisplay -nodesktop -nosplash -r "run_through_folders;exit"