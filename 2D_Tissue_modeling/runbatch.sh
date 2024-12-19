sbatch run_2d.slurm 0.1 CaMKII_db 500
# sbatch run_2d.slurm 0.1 Default 500
# sbatch run_2d.slurm 0.1 CaMKII_inhb 500
# sbatch run_2d.slurm 0 Default 500
# sbatch run_2d.slurm 0 CaMKII_db 500
# sbatch run_2d.slurm 0 CaMKII_inhb 500

sbatch run_2d_CaMKII_exclusion_workstation.slurm 0.1 CaMKII_db 500 No_CaMKII_GapG
sbatch run_2d_CaMKII_exclusion_workstation.slurm 0.1 CaMKII_db 500 No_CaMKII_NaV
sbatch run_2d_CaMKII_exclusion_workstation.slurm 0.1 CaMKII_db 500 No_CaMKII_K1