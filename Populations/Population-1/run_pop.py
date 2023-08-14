

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np






# def run(Para1, Para2, Para3, S2_amp):
def run(Para1):

	P=['./main_HAM_Signalling_cvode_new', str(BCL), '0', str(ISO),  str(CaMKII_inhb), str(CaMKII_db)]
	for i in range(len(Para1)):
		P.append(str(Para1[i]))

	call(P, stdout=f)


para_set1 = np.loadtxt('para.log')


f = open("AP.log.dat.1Hz", "w+")
# f2 = open("para_log.dat", "w+")

BCL = int(sys.argv[1])
ISO = float(sys.argv[2])
CaMKII_inhb = int(sys.argv[3])
CaMKII_db = int(sys.argv[4])


for i in range(600):
	index = str(i)
	print ('processing ID ' +  index)
	run(para_set1[i]);
	call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);
	call('mv Restart_ICs/ICs.bin Restart_ICs/ICs.bin.'+index, shell=True);

