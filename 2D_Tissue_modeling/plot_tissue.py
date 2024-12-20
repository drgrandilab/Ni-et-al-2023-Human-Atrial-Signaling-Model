
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

fonts = 22;
fonts_title = 22
leng_fontsize = 22
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


def plot_tissue(folder, ID):

	fig = plt.figure(figsize=(5.1,5))
	num_row = 1;
	num_col = 1;
	gs1 = gridspec.GridSpec(num_row,num_col);
	panel = {};
	
	for i in np.arange(num_row):
		for j in np.arange(num_col):
			panel[i,j] = plt.subplot(gs1[i,j]);
			ax = panel[i,j];
	
	
	NX = 125
	NY = 120
	
	index = ID
	# folder = 'reduce_diffusion/BCL.1000.ISO.0.1.CaMKII.CaMKII_db/'
	file = 'v_%04d.bin'%index
	data = np.fromfile(folder+file, dtype=np.float32);
	data = np.reshape(data, (NY, NX))
	
	
	ax = panel[0,0]
	
	ax.imshow(data, aspect='equal',cmap=plt.get_cmap('hot'), interpolation='gaussian',vmax=25, vmin=-88,origin='lower')#, vmax=0.4)
	ax.spines['left'].set_linewidth(0)
	ax.spines['right'].set_linewidth(0)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.yaxis.set_tick_params(width=0)
	ax.xaxis.set_tick_params(width=0)
	# ax.set_xlabel([])
	ax.set_xticklabels('')
	ax.set_yticklabels('')
	plt.tight_layout()
	plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99,
	wspace=0.28, hspace=0.10)
	plt.savefig(folder.replace('/', '_').replace('..', '')+str(index)+'.png', dpi=300)
	# plt.show()
	plt.close('all')




# for i in [2010,2040,2055,2350,2665]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2010,2040,2055,2350,2665]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2700]2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2010,2040,2055,2350,2665, 3200]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2010,2040,2055,2525,2610, 2625, 2663]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2662, 2675]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [3255, 3265, 3315, 3900]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
# for i in [2010,2040,2055,2525,2610, 2625, 2657]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:
for i in [3250]:#2010,2040,2055,2350,2665, 2680,2703,2705,2740,3350]:#[2010, 2040, 2055, 2070, 2100, 2200, 2300,2350, 2500, 2525, 2600, 2610, 2612, 2652, 2663, 2665, 2680, 3205, 2700, 2703,2705, 3200, 3203, 3220, 3255,3271, 3315, 3780, 3900]:

	print(i)
	# plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_inhb/', i)
	# plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_db/', i)
	# plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG/', i)
	# plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG/', i)
	plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV/', i)
	# plot_tissue('../BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1/', i)