

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np






def run(Para1):
	P=['./main_HAM_Signalling_cvode_phos_level', '500', '0', '0.1', '0', '1']
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	call(P, stdout=f)


def run_Restart(Para1, Para2, Para3, S2_amp):

	P=['./NCZ_Model_Threshold_long_5ms', '1000', '14', '1000', 'WT', 'Restart', '0', '0']
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	for i in range(len(Para2)):
		P.append(str(Para2[i]))
	for i in range(len(Para3)):  
		P.append(str(Para3[i]))
	P.append(str(S2_amp))
	call(P, stdout=f)





para_set1 = np.loadtxt('para_level.log')


f = open("AP.log.dat.1Hz", "w+")




for i in range(600):
	index = str(i)
	print ('processing ID ' +  index)
	run(para_set1[i]);
	call('mv HAM_wrap_out.dat AP.BCL.500.ID.'+index, shell=True);
	call('mv Restart_ICs/ICs.bin Restart_ICs/ICs.bin.'+index, shell=True);

