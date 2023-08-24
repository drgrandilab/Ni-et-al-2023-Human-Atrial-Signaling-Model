

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np






# def run(Para1, Para2, Para3, S2_amp):
def run(Para1):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./main_HAM_Signalling_cvode_new', str(BCL), '0', str(ISO),  str(CaMKII_inhb), str(CaMKII_db)]
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	# for i in range(len(Para2)):
	# 	P.append(str(Para2[i]))
	# for i in range(len(Para3)):  
	# 	P.append(str(Para3[i]))
	# P.append(str(-20))
	call(P, stdout=f)


para_set1 = np.loadtxt('para.log')


f = open("AP.log.dat.1Hz", "w+")
# f2 = open("para_log.dat", "w+")

BCL = int(sys.argv[1])
ISO = float(sys.argv[2])
CaMKII_inhb = int(sys.argv[3])
CaMKII_db = int(sys.argv[4])
# run(BCL,"Normal", 0,0,0,0);
# Mode="SimDrug"

# sys.stdout = open('file', 'w')

for i in range(600):
	index = str(i)
	print ('processing ID ' +  index)
	run(para_set1[i]);
	call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);
	call('mv Restart_ICs/ICs.bin Restart_ICs/ICs.bin.'+index, shell=True);

	# call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);
	# call('mv APD_measure.dat APD_measure.dat.ID.'+index, shell=True);
	# call('mv betaAR_out.dat betaAR_out.dat.ID.'+index, shell=True);
	# call('mv CaM_cyto_out.dat CaM_cyto_out.dat.ID.'+index, shell=True);
	# call('mv CaM_dyad_out.dat CaM_dyad_out.dat.ID.'+index, shell=True);
	# call('mv CaMKII_out.dat CaMKII_out.dat.ID.'+index, shell=True);
	# call('mv para_out.dat para_out.dat.ID.'+index, shell=True);