

import matplotlib

# matplotlib.use('Qt4Agg') 

import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
import csv
from matplotlib import rc
import pandas as pd

import numpy as np



def readlinedata(file):
	data = np.fromfile(file, sep = '',dtype=np.float32)

	# data = np.fromfile()

	# print(len(data))
	# print (data)

	nx = 55;
	ny = 17;
	nz = 11;
	# NZ = 
	linedata_longti = []
	linedata_trans = []
	for k in range(nz):
		for i in range(ny):
			for j in range(nx):
				if k == 5 and i == 9:
				# if j == 17 and i == 4:
					id = j + i * nx + k * nx * ny;
					linedata_longti.append(data[id])
				if k == 6 and j == 23:
				# if j == 17 and i == 4:
					id = j + i * nx + k * nx * ny;
					linedata_trans.append(data[id])
	# print(linedata_longti)
	return linedata_longti,linedata_trans

data = []
data_trans = []
# folder = 'Data_50%_Ncx'
# folder ='Data_Vp_x_2'
# folder ='Data'
# folder = 'Data_NCX_50%';start=14250
# folder = 'Data_NCX_50%_gamma_20%'
folder = 'AF.BCL.250.CKII_db.ISO.0';start=10 
# folder = 'AF.BCL.250.CKII_db.ISO.0';start=10 
folder = 'AF.CL.250.CKII_inh.ISO.0';start=10 
# folder = 'AF.CL.250.CKII_inh.ISO.1';start=10 
# folder = 'Data_CaT_50%'

# start=10820-500
# end = 11870+500`



end = start+15000
for i in range(start, end,10):
	file = folder+'/BinaryFiles/map_ci_%d.ID.0.bin'%i
	a,b=readlinedata(file)
	data.append([])
	data_trans.append([])
	data[-1]= a
	data_trans[-1]= b
# vm = np.loadtxt(folder+'/ci.txt',unpack=True)
## rcParams are the default parameters for matplotlib
import matplotlib as mpl
rc('mathtext',default='regular')
fonts = 16+8+8;
fonts_title = 16+8+8
leng_fontsize = 17+8
mpl.rcParams['font.size'] = fonts
# mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = fonts
mpl.rcParams['xtick.labelsize'] = fonts
mpl.rcParams['ytick.labelsize'] = fonts
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['ytick.major.pad']='8'
mpl.rcParams['xtick.major.pad']='10'



fig = plt.figure(figsize=(19*1.5/2.54,19*1/2.54))
num_row = 2;
num_col = 1;
gs1 = gridspec.GridSpec(num_row,num_col
	# width_ratios=[1.15,1],
    ,height_ratios=[1,1]
	);
panel = {};

for i in np.arange(num_row):
	for j in np.arange(num_col):
		if j>0:
			panel[i,j] = plt.subplot(gs1[i,j], sharey = panel[i,0]);
		else :
			panel[i,j] = plt.subplot(gs1[i,j]);

		ax = panel[i,j];

ax = panel[0,0]
# plt.figure(figsize=(20, 10))
# plt.subplot(411)
# ax.plot(vm[0]-start, vm[1],'k',lw=2*2)
# ax.set_xlim(start-start,end-start)
# ax.set_ylim(-70,20)
# ax.set_yticks([-40,0])
# ax.set_ylabel('V$_m$ (mV)')
ax.imshow(np.array(data).transpose(), aspect='auto',cmap=plt.get_cmap('hot'), interpolation='gaussian')#, vmax=0.4)#

ax = panel[1,0]

im=ax.imshow(np.array(data_trans).transpose(), aspect='auto',cmap=plt.get_cmap('hot'), interpolation='gaussian', vmax=0.9,vmin=0.15)#, )#

# plt.colorbar(im)

# ax = panel[2,0]


# # print(len(vm[0]))
# # ax.plot(np.mean(np.array(data), axis=1), 'k',lw=2*2)
# ax.set_xlim([0, len(data)])
# ax.set_ylabel('Ca$^{2+}$ ($\mu$M)')

# ax.set_ylim(0.1,0.4)
# ax.set_yticks([0.1,0.4])
# ax2 = ax.twinx()
# # plt.plot(np.mean(np.array(data), axis=1))
# # ax2.plot(vm[0]-start, vm[25-1],'0.5',alpha=0.8,lw=2*2)
# # ax2.set_xlim(start,end)
# ax2.set_ylabel('I$_{NCX}$ (pA/pF)')
# ax2.set_ylim(-3.5,0)
# ax2.set_yticks([-3,0])
# ax2.yaxis.set_tick_params(width=4)
# ax2.yaxis.set_tick_params(length=2*6)
# ax2.spines['right'].set_linewidth(8)
# ax2.spines['top'].set_linewidth(0)
# ax2.spines['bottom'].set_linewidth(0)
# ax.plot(np.mean(np.array(data), axis=1), 'k',lw=2*2)
# ax = panel[3,0]

# # plt.plot(np.mean(np.array(data), axis=1))
# ax.plot(vm[0]-start, vm[25-1],'k',lw=2*2)
# ax.set_xlim(start-start,end-start)
# # ax.xlim([0, len(data)])
# ax.set_ylabel('I$_{NCX}$ (pA/pF)')
# ax.set_xlabel('Time (ms)')
# ax.set_ylim(-3.5,0)
# ax.set_yticks([-3,0])

# for i in np.arange(num_col*num_row):
# 	ax = fig.axes[i]
# 	# ax.autoscale(enable=True, axis='x')  
# 	# ax.autoscale(enable=True, axis='r')  
# 	ax.spines['right'].set_visible(False)
# 	#ax.spines['left'].set_visible(False)
# 	ax.spines['top'].set_visible(False)
# 	ax.spines['bottom'].set_visible(True)
# 	ax.xaxis.set_ticks_position('bottom')
# 	ax.spines['bottom'].set_linewidth(8)
# 	ax.spines['left'].set_linewidth(8)
# 	ax.spines['right'].set_linewidth(8)
# 	ax.yaxis.tick_left()
# 	# if case1:
# 	# 	ax.set_xlim(24000,24900) 
# 	# else:
# 	# 	ax.set_xlim(22500, 24500) 

# 	# ax.set_ylim(-85,46)
# 	# # plt.xticks(np.linspace(0,600, 5))
# 	# ax.set_xticks([0, 250, 500])
# 	ax.yaxis.set_tick_params(width=4)
# 	ax.xaxis.set_tick_params(width=4)
# 	ax.yaxis.set_tick_params(length=2*6)
# 	ax.xaxis.set_tick_params(length=2*6)
# 	plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95,
# 		wspace=0.32, hspace=0.20)
# 	ax.yaxis.set_label_coords(-0.05, 0.5)
# for i in np.arange(3):
# 	ax = fig.axes[i]
# 	ax.set_xticklabels([])
# 	ax.spines['bottom'].set_visible(False)
# 	ax.xaxis.set_tick_params(width=0)

plt.savefig(folder+'.fig.pdf')
plt.show()




