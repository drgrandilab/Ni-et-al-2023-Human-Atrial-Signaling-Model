

run_single () {

BCL=$1
ISO=$2
CaMKII_inhb=$3
CaMKII_db=$4

folder=BCL.$1.ISO.$2.CaMKII_inhb.$3.CaMKII_db.$4
mkdir -p $folder


cp para_large_sigma.log.20000 run.slurm run_pop.py main_HAM_Signalling_cvode_new $folder
cp -r Restart_ICs $folder
cd $folder 
sbatch run.slurm $BCL $ISO $CaMKII_inhb $CaMKII_db
cd ../
}

# make clean;
# make


# run_single 1000 0.1 0 1
# run_single 1000 0.1 1 0
# run_single 1000 0.1 0 0
# run_single 1000 0 0 0
# run_single 1000 0 1 0
# run_single 1000 0 0 1

run_single 500 0.1 0 1
run_single 500 0.1 0 0
run_single 500 0.1 1 0
run_single 500 0 1 0
run_single 500 0 0 1
run_single 500 0 0 0



# run_single 1000 0 0 0
