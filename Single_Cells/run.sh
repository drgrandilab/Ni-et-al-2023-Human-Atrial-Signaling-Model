#bash

run_stuff() {
	./main_HAM_Signalling_cvode_new $1 $2 $3 $4 $5

	DAT=(APD_measure.dat HAM_wrap_out.dat para_out.dat) # betaAR_matlab.dat betaAR_out.dat  CaM_cyto_out.dat  CaM_dyad_out.dat  CaMKII_out.dat  CaM_sl_out.dat    LTCC_current.dat  para_out.dat)

	for ((i = 0; i < ${#DAT[@]}; ++i));
	do
		mv ${DAT[$i]} ${DAT[$i]}.$1.AF.$2.ISO.$3.CKII_ih.$4.CKII_db.$5
	done

 # betaAR_out.dat  CaM_cyto_out.dat  CaM_dyad_out.dat  CaMKII_out.dat  CaM_sl_out.dat  HAM_wrap_out.dat  LTCC_current.dat  para_out.dat

}

run_stuff 1000 0 0 0 0 >> APD_out.dat
run_stuff 1000 0 0.1 0 0 >> APD_out.dat
run_stuff 333.3333 0 0 0 0 >> APD_out.dat


