#!/bin/bash
#SBATCH --job-name=1Hz
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# bash run.sh
python run_pop.py $1 $2 $3 $4
