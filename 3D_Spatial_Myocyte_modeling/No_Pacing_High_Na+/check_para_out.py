
import numpy as np


import matplotlib.pyplot as plt

from math import *

def check_paraout(file):

	data = np.loadtxt(file, unpack=True)

	print(file)

	dct = {}


	for a,b in zip(['LCC_CKdyadp', 'RyR_CKp', 'PLB_CKp', 'NaV_CKp', 'Ikur_CKp', 'Itof_CKp'], range(6)):
		print(f"{a} { data[b+1].mean():.4f}")

		dct[a]=data[b+1].mean()

	dct['Gapjunct_CKp'] = dct['NaV_CKp'] 
	dct['IK1_CKp'] = dct['NaV_CKp'] 
	# pass

	Nav_CKp_wt = 2.5654 / 30; # 8.5# 1X, 1 Hz
	Nav_CKp_oe = 23.9069 / 30; # 80# 6X, 1 Hz

	Ikur_CKp_ko = 0;
	Ikur_CKp_wt = 2.5654 / 30;

	Itof_CKp_wt = 0.11;#//2.5654 / 30.0;

	fckiim2_j = dct['LCC_CKdyadp'] * 0.1
	fCKII_RyR = 1 + 3 * (dct['RyR_CKp'] - 0.21);
	fCKII_PLB = (1 - 0.5 * dct['PLB_CKp'] );

	kCKII_Nav_change = (dct['NaV_CKp'] - Nav_CKp_wt) / (Nav_CKp_oe - Nav_CKp_wt); # 0 WT, 1 OE
	kCKII_Ikur = 0.5 + 1.0 / ( 1 + exp( (dct['Ikur_CKp'] - Ikur_CKp_wt) / -0.08 ));


	kCKII_Itof_vshift = -5 + 10 / ( 1 + exp( (dct['Itof_CKp'] - Itof_CKp_wt) / -0.1 ));
	kCKII_Itof_tau = 0.5 + 1.0 / ( 1 + exp( (dct['Itof_CKp'] - Itof_CKp_wt) / -0.1 )) ;
	kCKII_IK1_G = 0.5 + 1.0 / ( 1 + exp( (dct['IK1_CKp']  - 0.1) / -0.15 ))

	kCKII_Gapjunct_G = 1.35 - (0.9 / (1 + exp( -(dct['Gapjunct_CKp'] - 0.16) / 0.12 ) ) );
	fPKA_RyR = 1.0
	RyR_koSRCa_Scale = fCKII_RyR + fPKA_RyR - 1.0;



files = ['Data_BCL.250.ISO.0.CaMKII_inh.0.CaMKII_db.0/para_out.dat',
'Data_BCL.333.333.ISO.0.CaMKII_inh.0.CaMKII_db.0/para_out.dat',
'Data_BCL.500.ISO.0.CaMKII_inh.0.CaMKII_db.0/para_out.dat',
'Data_BCL.1000.ISO.0.CaMKII_inh.0.CaMKII_db.0/para_out.dat',
'Data_BCL.2000.ISO.0.CaMKII_inh.0.CaMKII_db.0/para_out.dat']

files_db = ['Data_BCL.250.ISO.0.CaMKII_inh.0.CaMKII_db.1/para_out.dat',
'Data_BCL.333.333.ISO.0.CaMKII_inh.0.CaMKII_db.1/para_out.dat',
'Data_BCL.500.ISO.0.CaMKII_inh.0.CaMKII_db.1/para_out.dat',
'Data_BCL.1000.ISO.0.CaMKII_inh.0.CaMKII_db.1/para_out.dat',
'Data_BCL.2000.ISO.0.CaMKII_inh.0.CaMKII_db.1/para_out.dat']


for file in files_db:
	check_paraout(file)
