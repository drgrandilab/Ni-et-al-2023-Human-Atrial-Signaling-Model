

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np

import os
from joblib import Parallel, delayed
import multiprocessing




# def run(Para1, Para2, Para3, S2_amp):
def run(Para1, Popul_ID):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./main_HAM_Signalling_cvode_new', str(BCL), str(Popul_ID), str(ISO),  str(CaMKII_inhb), str(CaMKII_db)]
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	call(P, stdout=f)


Popul_params = np.loadtxt('para_large_sigma.log.20000')



# f2 = open("para_log.dat", "w+")

BCL = float(sys.argv[1])
ISO = float(sys.argv[2])
CaMKII_inhb = int(sys.argv[3])
CaMKII_db = int(sys.argv[4])
# run(BCL,"Normal", 0,0,0,0);
# Mode="SimDrug"
f = open("AP.log.dat.20000", "w+")
# sys.stdout = open('file', 'w')

IDs_to_run = list(range(len(Popul_params)))


run_parallel = True #False

if not run_parallel:
	for i in IDs_to_run: #range(600):
		index = str(i)
		print ('processing ID ' +  index)
		run(Popul_params[i], i);
		call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);
		call('mv Restart_ICs/ICs.bin Restart_ICs/ICs.bin.'+index, shell=True);



if run_parallel:
	num_cores = 40#multiprocessing.cpu_count() - 2

	results=Parallel(n_jobs=num_cores, prefer="threads")(
	 		delayed(run) 
	 				(Popul_params[PopulID], PopulID) 
	 				for PopulID in IDs_to_run)

f.close()