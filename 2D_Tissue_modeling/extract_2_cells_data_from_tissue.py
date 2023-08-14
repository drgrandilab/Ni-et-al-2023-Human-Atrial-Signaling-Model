import numpy as np


def readcelldata(file, index):
	# data = np.fromfile(file, sep = '\t')
	tmp = np.fromfile(file, dtype=np.float32) 
	
	# data = np.fromfile()
	
	# print(len(data))
	# print (data)
	
	# nx = 107;
	# ny = 101;
	# # nz = 11;
	# # NZ = 
	# linedata = []

	# for i in range(ny):
	# 	for j in range(nx):
	# 		if k == 5 and j == 25:
	# 			id = j + i * nx + k * nx * ny;
	# 			linedata.append(data[id])
	# print(linedata)
	linedata = tmp[index]
	return linedata







Cai = []
Cai2 = []
NX=125

ID1=20*NX + 10
ID2=20*NX + 80
# # ID2 =5000+350+500+9
# # folder = 'BCL.500.ISO.0.1.CaMKII.CaMKII_inhb'
# # folder = 'BCL.500.ISO.0.CaMKII.CaMKII_db'
# folder = 'BCL.1000.ISO.0.1.CaMKII.CaMKII_db'
# # folder = 'BCL.500.ISO.0.1.CaMKII.Default'
# # folder = 'BCL.500.ISO.0.CaMKII.Default'
# # folder = 'BCL.1000.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG'
# folder = 'BCL.1000.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV'
# # folder = 'BCL.1000.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1'


def output_2Cells(folder):
	data = []
	data2 = []
	for i in range( 1,10000,1):
		file = folder+'/v_%04d.bin'%i
		data.append([])
		data2.append([])
		data[-1]= readcelldata(file,ID1)
		data2[-1]= readcelldata(file,ID2)
		# Cai.append([])
		# Cai2.append([])

		# file = 'BINFILEs/Cai_%04d.bin'%i
		# Cai[-1]= readcelldata(file,ID1)
		# Cai2[-1]= readcelldata(file,ID2)



	np.savetxt(folder+'.data1.dat', data, fmt='%.5f')
	np.savetxt(folder+'.data2.dat', data2, fmt='%.5f')


folders = [
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG',
'BCL.500.ISO.0.CaMKII.CaMKII_inhb',
'BCL.500.ISO.0.CaMKII.CaMKII_db',
'BCL.500.ISO.0.CaMKII.Default',
'BCL.500.ISO.0.1.CaMKII.CaMKII_inhb',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db',
'BCL.500.ISO.0.1.CaMKII.Default'
]


foldersa = [
# 'BCL.500.ISO.0.CaMKII.CaMKII_inhb',
# 'BCL.500.ISO.0.CaMKII.CaMKII_db',
# 'BCL.500.ISO.0.CaMKII.Default',
# 'BCL.500.ISO.0.1.CaMKII.CaMKII_inhb',
# 'BCL.500.ISO.0.1.CaMKII.CaMKII_db',
# 'BCL.500.ISO.0.1.CaMKII.Default',
# 'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG',
# 'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1',

]



for i in folders:
	print(i)
	output_2Cells(i)

# import numpy as np
# import collections
# import matplotlib.gridspec as gridspec
# import matplotlib.pyplot as plt
# import csv
# from numpy import *

# from matplotlib import rc
# import pandas as pd
# ### rcParams are the default parameters for matplotlib
# import matplotlib as mpl
# rc('mathtext',default='regular')

# fonts = 17;
# fonts_title = 17
# leng_fontsize = 17
# mpl.rcParams['font.size'] = fonts
# mpl.rcParams['font.family'] = 'Arial'
# mpl.rcParams['axes.labelsize'] = fonts
# mpl.rcParams['xtick.labelsize'] = fonts
# mpl.rcParams['ytick.labelsize'] = fonts
# mpl.rcParams['axes.linewidth'] = 1.
# mpl.rcParams['ytick.major.pad']='8'
# mpl.rcParams['xtick.major.pad']='10'

# Q10 = 2.0

# Q10_Vact_Shift = 4.0
# Q10_Inaact_Shift = 4.0
# Temp = 22.0
# # rc('text', usetex=True)
# # rc('font',**{'family':'Times','sans-serif':['aria']})

# xpos = -0.4
# ypos = 0.5
# lblsize = 22
# boldness = 500
# ticklblsize = 18
# transparency = 1.0

# fig = plt.figure(figsize=(8.4*1.1,5.5))
# num_row = 1;
# num_col = 1;
# gs1 = gridspec.GridSpec(num_row,num_col);
# panel = {};

# for i in np.arange(num_row):
# 	for j in np.arange(num_col):
# 		panel[i,j] = plt.subplot(gs1[i,j]);
# 		ax = panel[i,j];


# color = ['k','r']




# col = 0;
# row = 0;
# ax = panel[row,col];
# #ax.plot(data0[0]+6e4, data0[1], 'k',lw= 2.0, label='Original')
# # ax.plot(data4[0], data4[1], 'r--',lw= 2.0, label='PCL = 250')
# ax.plot(data2, color[1],lw= 3, alpha=1, label='Baseline')
# ax.plot(data, color[0],lw= 3, alpha=1, label='Baseline')

# ax.set_ylabel('Voltage (mV)')
# #ax.set_xlabel('Time (ms)')
# # ax.set_yticks([-80,-40,0,20])
# # ax.set_xlim([0,3000])

# col = 0;
# row = 1;
# # ax = panel[row,col];

# # #ax.plot(data0[0]+6e4, data0[7], color[i],lw= 1.5, alpha=1, label='Original')
# # # ax.plot(Cai, color[0],lw= 1.5, alpha=1, label='Baseline')
# # # ax.plot(Cai2, color[1],lw= 1.5, alpha=1, label='Baseline')
# # ax.plot(np.mean(np.array(data_cai_5000), axis=1), lw=3, c=color[1])
# # ax.plot(np.mean(np.array(data_cai_5), axis=1), lw=3, c=color[0])

# # ax.set_ylabel('$[Ca^{2+}]_i$ ($\mu$M)')
# # ax.set_yticks([0.1,0.3,0.5])
# # ax.set_xlim([0,3000])


# # ax = panel[2, 0]
# # ax.imshow(np.array(data_cai_5000).transpose(), aspect='auto',cmap=plt.get_cmap('hot'), interpolation='none', vmax=0.8)
# # ax.spines['left'].set_linewidth(0)
# # ax.spines['right'].set_linewidth(0)
# # ax.spines['right'].set_visible(False)
# # ax.spines['left'].set_visible(False)
# # ax.spines['top'].set_visible(False)
# # ax.spines['bottom'].set_visible(False)
# # ax.yaxis.set_tick_params(width=0)
# # ax.xaxis.set_tick_params(width=0)

# # ax = panel[3, 0]
# # ax.imshow(np.array(data_cai_5).transpose(), aspect='auto',cmap=plt.get_cmap('hot'), interpolation='none', vmax=0.8)
# # ax.spines['left'].set_linewidth(0)
# # ax.spines['right'].set_linewidth(0)
# # ax.spines['right'].set_visible(False)
# # ax.spines['left'].set_visible(False)
# # ax.spines['top'].set_visible(False)
# # ax.spines['bottom'].set_visible(False)
# # ax.yaxis.set_tick_params(width=0)
# # ax.xaxis.set_tick_params(width=0)


# for i in range(1):
# 	ax = fig.axes[i]
# 	# ax.spines['left'].set_position(('data', -100))
# 	ax.spines['right'].set_color('none')
# 	# ax.spines['bottom'].set_position(('data', 0))
# 	ax.spines['top'].set_color('none')
# 	# ax.spines['left'].set_smart_bounds(True)
# 	# ax.spines['bottom'].set_smart_bounds(True)
# 	ax.xaxis.set_ticks_position('bottom')
# 	ax.yaxis.set_ticks_position('left')
# 	ax.xaxis.set_tick_params(width=2)
# 	ax.yaxis.set_tick_params(width=2)

# 	# for i in range(0, len(sys.argv)-1):
# 	ax.spines['bottom'].set_linewidth(2)
# 	ax.spines['left'].set_linewidth(2)
# 	ax.spines['right'].set_linewidth(2)
# 	# 	norm = abs(min(IV[i][1]))
# 	# 	ax.plot(IV[i][0], IV[i][1]/norm, label=labels[i], lw=2)
# 	# ax.set_xlim(4000,20000)
# 	# ax.set_xlim(101200,113325)
# 	# ax.set_xticks(np.linspace(1,20,6))
# 	plt.subplots_adjust(left=0.12, bottom=0.02, right=0.99, top=0.99,
# 		wspace=0.28, hspace=0.10)
# 	ax.yaxis.set_label_coords(-0.08, 0.5)


# for i in range(1):
# 	ax = fig.axes[i]
# 	# ax.spines['left'].set_position(('data', -100))
# 	ax.spines['right'].set_color('none')
# 	# ax.spines['bottom'].set_position(('data', 0))
# 	ax.spines['bottom'].set_color('none')
# 	ax.set_xticks([])
# 	ax.xaxis.set_tick_params(width=0)
# 	# ax.yaxis.set_tick_params(width=0)
# 	# ax.spines['left'].set_smart_bounds(True)
# # plt.savefig('Ito_react.pdf')
# plt.savefig('AP_Cai.pdf')

# plt.show()

