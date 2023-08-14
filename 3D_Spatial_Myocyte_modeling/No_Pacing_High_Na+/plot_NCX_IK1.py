import numpy as np
import collections
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import csv
from numpy import *

from matplotlib import rc
import pandas as pd
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




# folder = 'BCL.500.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1'
# folder = 'BCL.500.ISO.0.1.CaMKII_inhb.0.CaMKII_db.0'
# folder = 'BCL.500.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0'

folder = 'BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.1'
folder = 'BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.0'
folder = './'
files = ['AF.BCL.250.CKII_db.ISO.0/global_result/global_current.txt','AF.CL.250.CKII_inh.ISO.0/global_result/global_current.txt']

fig = plt.figure(figsize=(8.4/1.8,6*1.1/1.2))
num_row = 3;
num_col = 1;
gs1 = gridspec.GridSpec(num_row,num_col);
panel = {};


colors = ['tab:red','tab:blue'];


index = [1,2]
for ci,f,index_i in zip(colors,files,index):


	for i in np.arange(num_row):
		for j in np.arange(num_col):
			panel[i,j] = plt.subplot(gs1[i,j]);
			ax = panel[i,j];

	# data0 = np.loadtxt('Grandi_wrap_out.dat.OriModel', unpack=True)
	# data1 = np.loadtxt('CM_0_0/AP.BCL.1000.ID.%d'%i, unpack=True)
	data1=[]
	data1 = np.array(pd.read_csv(f, sep = '\t')).transpose()
	# data2 = np.loadtxt('CM_0_0/', unpack=True)
	# data3 = np.loadtxt('with_CMKIII/new/HAM_wrap_out.dat.1000.AF.0.ISO.0.1.IKur.1.CKII_ih.1.ISK.0.5', unpack=True)
	# data4 = np.loadtxt('ISO_CKI/HAM_wrap_out.dat.1000.AF.0.ISO.0.1.IKur.1.CKII_ih.1', unpack=True)
	
	col = 0;
	row = 0;
	ax = panel[row,col];
	#ax.plot(data0[0]+6e4, data0[1], 'k',lw= 2.0, label='Original')
	# ax.plot(data4[0], data4[1], 'r--',lw= 2.0, label='PCL = 250')
	ax.plot(data1[0], data1[1], ci,lw= 1.5, alpha=1, label='Baseline')
	# ax.plot(data3[0], data3[1], 'b',lw= 1.5, alpha=0.2, label='PCL = 300')
	# ax.plot(data2[0], data2[1], 'r',lw= 1.5, alpha=0.2, label='ISO')
	# ax.set_xscale("log")
	
	# V=-80
	# print Ito_kv43_inActime(V), Ito_kv14_inactimss(V), Ito_kv14_inactim(V), ICaL_f_tau(V)
	
	# ax.legend()
	# leg = ax.legend(numpoints = 1, bbox_to_anchor=(0.8, 0.8), fancybox=False, prop={'size':16}, labelspacing=0.3,handletextpad=1, handlelength=2)
	# leg.draggable(state=True)
	# leg.draw_frame(False)
	# leg.get_frame().set_alpha(0.0)
	ax.set_ylabel('V$_m$ (mV)')
	#ax.set_xlabel('Time (ms)')
	ax.set_yticks([-80,-40,0,40])
	# ax.set_ylim([-90., 40])

	
	col = 0;
	row = 1
	if 	row == index_i:
		ax = panel[row,col];
		
		#ax.plot(data0[0]+6e4, data0[7], 'k',lw= 1.5, alpha=0.2, label='Original')
		ax.plot(data1[0], data1[20], "tab:red",lw= 1.5, alpha=1, label='PCL = 1000')
		ax.plot(data1[0], data1[5], "tab:red",lw= 1.5, alpha=1, label='PCL = 1000')
		# ax.plot(data3[0], 1000*data3[7], 'b',lw= 2.0, label='PCL = 300')
		# ax.plot(data2[0], 1000*data2[7], 'r',lw= 2.0, label='PCL = 500')
		# ax.plot(data4[0], 1000*data4[7], 'r--',lw= 2.0, label='PCL = 250')
		# ax.set_ylabel('$[Ca^{2+}]_i$ ($\mu$M)')
		# ax.set_yticks(np.linspace(0.,0.4,3))
		# ax.plot([600500-50+100, 600500-50 +100+1000], [0.05, 0.05], lw =3, c=ci)
		ax.set_ylim([-3, 0.4 ])
		ax.set_yticks(np.linspace(-3,0,3))
		# ax.plot([310e3-100+200,310e3-100+200+2e3], [1,1], lw=4, color=ci)
		ax.set_ylabel('(pA/pF)')
	

	col = 0;
	row = 2
	if 	row == index_i:
		ax = panel[row,col];
		
		
		#ax.plot(data0[0]+6e4, data0[7], 'k',lw= 1.5, alpha=0.2, label='Original')
		ax.plot(data1[0], data1[20], "tab:blue",lw= 1.5, alpha=1, label='PCL = 1000')
		ax.plot(data1[0], data1[5], "tab:blue",lw= 1.5, alpha=1, label='PCL = 1000')
		# ax.plot(data3[0], 1000*data3[7], 'b',lw= 2.0, label='PCL = 300')
		# ax.plot(data2[0], 1000*data2[7], 'r',lw= 2.0, label='PCL = 500')
		# ax.plot(data4[0], 1000*data4[7], 'r--',lw= 2.0, label='PCL = 250')
		# ax.set_ylabel('$[Ca^{2+}]_i$ ($\mu$M)')
		ax.set_yticks(np.linspace(0.,0.4,3))
		# ax.plot([600500-50+100, 600500-50 +100+1000], [0.05, 0.05], lw =3, c=ci)
		ax.set_ylim([-3, 0.4 ])
		ax.set_yticks(np.linspace(-3,0,3))
		# ax.plot([310e3-100+200,310e3-100+200+2e3], [1,1], lw=4, color=ci)
		ax.set_ylabel('(pA/pF)')
		# ax.plot(data3[0], 1000*data3[7], 'b',lw= 2.0, label='PCL = 300')
		ax.plot([3000,4000], [-2,-2], lw =3, c='k')
		ax.set_ylim([-3, 0.4 ])
		ax.set_yticks(np.linspace(-3,0,3))


	for i in range(3):
		ax = fig.axes[i]
		# ax.spines['left'].set_position(('data', -100))
		ax.spines['right'].set_color('none')
		# ax.spines['bottom'].set_position(('data', 0))
		ax.spines['top'].set_color('none')
		# ax.spines['left'].set_smart_bounds(True)
		# ax.spines['bottom'].set_smart_bounds(True)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_tick_params(width=2)
		ax.yaxis.set_tick_params(width=2)

		# for i in range(0, len(sys.argv)-1):
		ax.spines['bottom'].set_linewidth(2)
		ax.spines['left'].set_linewidth(2)
		ax.spines['right'].set_linewidth(2)
		# 	norm = abs(min(IV[i][1]))
		# 	ax.plot(IV[i][0], IV[i][1]/norm, label=labels[i], lw=2)
		ax.set_xlim(-200,5000)
		# ax.set_xticks(np.linspace(1,20,6))
		plt.subplots_adjust(left=0.18, bottom=0.05, right=0.99, top=0.95,
			wspace=0.28, hspace=0.25)
		ax.yaxis.set_label_coords(-0.15, 0.5)


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



plt.savefig('_test_noISO_IK1_INCX.pdf', dpi=300)

plt.show()


