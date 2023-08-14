import numpy as np
import collections
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import csv
from numpy import *

from matplotlib import rc
# from load import *

### rcParams are the default parameters for matplotlib
import matplotlib as mpl
rc('mathtext',default='regular')

fonts = 14;
fonts_title = 14
leng_fontsize = 14
mpl.rcParams['font.size'] = fonts
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = fonts
mpl.rcParams['xtick.labelsize'] = fonts
mpl.rcParams['ytick.labelsize'] = fonts
mpl.rcParams['axes.linewidth'] = 1.
mpl.rcParams['ytick.major.pad']='8'
mpl.rcParams['xtick.major.pad']='10'

Q10 = 2.0

Q10_Vact_Shift = 4.0
Q10_Inaact_Shift = 4.0
Temp = 22.0
# rc('text', usetex=True)
# rc('font',**{'family':'Times','sans-serif':['aria']})

xpos = -0.4
ypos = 0.5
lblsize = 22
boldness = 500
ticklblsize = 18
transparency = 1.0

fig = plt.figure(figsize=(8.4/1.2*1.3,6*1.33*1.2))
num_row = 6;
num_col = 1;
gs1 = gridspec.GridSpec(num_row,num_col);
panel = {};

for i in np.arange(num_row):
	for j in np.arange(num_col):
		panel[i,j] = plt.subplot(gs1[i,j]);
		ax = panel[i,j];

# data0 = np.loadtxt('Grandi_wrap_out.dat.OriModel', unpack=True)
data = []


folder = '../D_25%/'

files = [
'BCL.500.ISO.0.1.CaMKII.Default',
'BCL.500.ISO.0.1.CaMKII.CaMKII_inhb',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV',
'BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG',
]
# 'BCL.500.ISO.0.CaMKII.CaMKII_inhb',
# 'BCL.500.ISO.0.CaMKII.CaMKII_db',
# 'BCL.500.ISO.0.CaMKII.Default',
for file in files:
	data.append(np.loadtxt(folder+file+'60_0_0_.dat', unpack=True) )


# data.append ( np.loadtxt('BCL.500.ISO.0.CaMKII.Default.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.CaMKII.CaMKII_db.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.CaMKII.CaMKII_inhb.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.CaMKII.Default.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.1.CaMKII.Default.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.1.CaMKII.CaMKII_inhb.data2.dat', unpack=True) )
# data.append ( np.loadtxt('BCL.500.ISO.0.1.CaMKII.CaMKII_db.data2.dat', unpack=True) )
lables = files

# data2 = np.loadtxt('BCL.1000.ISO.0.CaMKII.CaMKII_db/data2.dat', unpack=True)
# data3 = np.loadtxt('BCL.1000.ISO.0.1.CaMKII.Default/data2.dat', unpack=True)
# data4 = np.loadtxt('BCL.1000.ISO.0.1.CaMKII.CaMKII_db/data2.dat', unpack=True)

for i in range(len(data)):
	ax = fig.axes[i]

	ax.plot(data[i][0],data[i][1], 'k',lw= 2.5, label=lables[i])
	ax.legend( frameon=False)


# # data_grandi = np.loadtxt('Grandi_wrap_out.dat', unpack=True);
# col = 0;
# row = 0;
# ax = panel[row,col];
# #ax.plot(data0[0]+6e4, data0[1], 'k',lw= 2.5, label='Original')
# # data_grandi
# # ax.plot(data_grandi[0]+450e3, data_grandi[1], 'grey',lw= 2.5,ls = '--', label='Grandi et al. 2011 Model')

# ax.plot(data1, 'k',lw= 2.5, label='PCL = 500')

# col = 0;
# row = 1;
# ax = panel[row,col];
# ax.plot(data2, 'k',lw= 2.5, label='PCL = 500')
# ax = panel[2,col];
# ax.plot(data3, 'k',lw= 2.5, label='PCL = 500')

# ax = panel[3,col];
# ax.plot(data4, 'k',lw= 2.5, label='PCL = 500')

for i in range(len(data)):
	ax = fig.axes[i]
	# ax.spines['left'].set_position(('data', -100))
	ax.spines['right'].set_color('none')
	# ax.spines['bottom'].set_position(('data', 0))
	ax.spines['top'].set_color('none')
	# ax.spines['left'].set_smart_bounds(True)
	# ax.spines['bottom'].set_smart_bounds(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_tick_params(width=2,length=4)
	ax.yaxis.set_tick_params(width=2,length=4)
	# ax.set_ylim(-90,15)

	# for i in range(0, len(sys.argv)-1):
	ax.spines['bottom'].set_linewidth(2)
	ax.spines['left'].set_linewidth(2)
	ax.spines['right'].set_linewidth(2)
	# 	norm = abs(min(IV[i][1]))
	# 	ax.plot(IV[i][0], IV[i][1]/norm, label=labels[i], lw=2)
	ax.set_xlim(-200,10000)
	# ax.set_xticks(np.linspace(1,20,6))
	plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95,
		wspace=0.28, hspace=0.30)
	ax.yaxis.set_label_coords(-0.1, 0.5)


for i in range(3):
	ax = fig.axes[i]
	# ax.spines['left'].set_position(('data', -100))
	ax.spines['right'].set_color('none')
	# ax.spines['bottom'].set_position(('data', 0))
	ax.spines['bottom'].set_color('none')
	ax.set_xticks([])
	ax.xaxis.set_tick_params(width=0)
	# ax.yaxis.set_tick_params(width=0)
	# ax.spines['left'].set_smart_bounds(True)
# plt.savefig('Ito_react.pdf')
plt.savefig('AP_ISO.pdf', dpi=300)

plt.show()



# fig = plt.figure(figsize=(8.4/1.5/5,6*1.5/5))
# num_row = 1;
# num_col = 1;
# gs1 = gridspec.GridSpec(num_row,num_col);
# panel = {};

# for i in np.arange(num_row):
# 	for j in np.arange(num_col):
# 		panel[i,j] = plt.subplot(gs1[i,j]);
# 		ax = panel[i,j];
# # ax.plot(data_grandi[0]+450e3, data_grandi[3], 'grey',lw= 2.,ls = '--', label='Grandi et al. 2011 Model')


# ax.plot(data2[0], data2[8], 'b',lw= 2.5, label='PCL = 500')
# ax.plot(data3[0], data3[8], 'g',lw= 2.5, label='PCL = 300')
# ax.plot(data4[0], data4[8], 'r',lw= 2.5, label='PCL = 250')
# ax.plot(data1[0], data1[8], 'k',lw= 2.5, label='PCL = 1000')
# ax.set_ylabel('I$_{CaL}$ (A/F)')
# # ax.set_yticks(np.linspace(-8,0,3))
# ax.set_ylim([-8,-5])
# ax.set_xlim(600000+4-1+2, 600000+10+5)
# plt.savefig('AP_compStep_insert.pdf', dpi=300)

# plt.show()
