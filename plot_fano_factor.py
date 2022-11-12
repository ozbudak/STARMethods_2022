"""
Plot mean expression vs. Fano factor 
Copyright (C) 2016 Ahmet Ay, Dong Mai, Soo Bin Kwon

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys, shared, os
import numpy, math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import xlrd, xlwt
from xlrd import XLRDError
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2,5,10,15,20]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 0.1

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=5)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def bin_fixedbinsize(mean, ff): # bin fano_factor into five groups based on mean (bin size is fixed)
	binned_mean = []
	binned_ff = []	
	for i in range(5):
		binned_mean.append([])
		binned_ff.append([])
		
	bin_size = (max(mean)-min(mean))/5 # same bin size for all bins
	for i in range(len(mean)):
		for j in range(5):
			if mean[i]<=min(mean)+bin_size*(j+1):
				binned_mean[j].append(mean[i])
				binned_ff[j].append(ff[i])				
				break	
				
	return binned_mean, binned_ff

def bin_fixednumslice(mean, ff): # bin fano_facot into five groups based on mean (number of slices in each bin is fixed)
	binned_mean = []
	binned_ff = []	
	for i in range(5):
		binned_mean.append([])
		binned_ff.append([])
		
	sorted_mean_her1 = sorted(mean) # rank the slices
	sorted_ff_her1 = [y for (x,y) in sorted(zip(mean,ff))]
	bin_size = int(math.floor(len(mean)/5))	# each bin gets the same number of slices (except for the last one)
	for i in range(4):
		binned_mean[i] = sorted_mean_her1[i*bin_size:(i+1)*bin_size]
		binned_ff[i] = sorted_ff_her1[i*bin_size:(i+1)*bin_size]	
	binned_mean[4] = sorted_mean_her1[4*bin_size:] # last bin takes the remaining data
	binned_ff[4] = sorted_ff_her1[4*bin_size:]	
	
	return binned_mean, binned_ff
	
def plot(ax, binned_mean, binned_ff, color): # plot mean expression vs. Fano factor
	xmax = float('-inf')
	ymax = float('-inf')		
	for i in range(5):
		p = ax.scatter(numpy.mean(binned_mean[i]),numpy.mean(binned_ff[i]), c=color, s = 30, edgecolors='none')
		ax.errorbar(numpy.mean(binned_mean[i]),numpy.mean(binned_ff[i]),
			xerr=2*numpy.std(binned_mean[i])/math.sqrt(len(binned_mean[i])),			
			yerr=2*numpy.std(binned_ff[i])/math.sqrt(len(binned_ff[i])), ls='none', c=color, capsize=4)			
		xmax = max(xmax, numpy.mean(binned_mean[i])+2*numpy.std(binned_mean[i])/math.sqrt(len(binned_mean[i])))
		ymax = max(ymax, numpy.mean(binned_ff[i])+2*numpy.std(binned_ff[i])/math.sqrt(len(binned_ff[i])))
	return p,xmax,ymax
			
def write(ws, binned_mean, binned_ff): # write data to Excel 	
	labels = ['Bin #','# slices','Mean RNA level','Std error','Fano factor','Std error']
	for i in range(len(labels)):
		ws.write(0,i,labels[i])	
	for i in range(5):
		if len(binned_mean[i])>0:
			line = [i+1,len(binned_mean[i]),numpy.mean(binned_mean[i]),numpy.std(binned_mean[i])/math.sqrt(len(binned_mean[i])),
				numpy.mean(binned_ff[i]),numpy.std(binned_ff[i])/math.sqrt(len(binned_ff[i]))]
		else:
			line = ["","","","",""]
		for j in range(len(line)):
			ws.write(i+1,j,line[j])
			
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('plot_fano_factor.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print ('plot_fano_factor.py: Number of embryos must be larger than zero.')
		exit(1)		
	
	num_embryos = int(sys.argv[1])
	    	        
	if len(sys.argv)==num_embryos+3:
		inputs = sys.argv[2:num_embryos+2]
		directory = sys.argv[num_embryos+2]		
	else:
		usage()			
		        	
	slice_mean_her1 = []
	slice_mean_her7 = []
	slice_fano_her1 = []
	slice_fano_her7 = []
          
        # Read and process input      
	for i in range(len(inputs)): # embryo
		# Open embryo data	
		if not os.path.isfile(inputs[i]):
			print ('plot_fano_factor.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('plot_fano_factor.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		worksheets = workbook.sheets()
		
		for j in range(len(worksheets)): # region 
			worksheet = worksheets[j]
			file_len =  worksheet.nrows	
								
			for k in range(1,file_len): # slice				
				row = list(worksheet.row(k))
					
				if isinstance(row[1].value,float): # valid slice
					num_cells = int(row[1].value) # number of cells in this slice
						
					her1 = []
					her7 = []
					for l in range(num_cells):
						# Take the cell's data only if its expression levels are positive after background subtraction	
						if row[8+2*l].value>0 and row[8+2*l+1].value>0:
							her1.append(row[8+2*l].value)
							her7.append(row[8+2*l+1].value)
					
					if len(her1)>=3: # valid only if 3 or more cells
						her1_mean = numpy.mean(her1)
						her7_mean = numpy.mean(her7)
																		
						slice_mean_her1.append(her1_mean) # her1 mean
						slice_mean_her7.append(her7_mean) # her7 mean
						
						# Intrinsic noise
						intrinsic_noise = 0						
						for m in range(len(her1)):
							intrinsic_noise += (her1[m]/her1_mean - her7[m]/her7_mean)**2
						intrinsic_noise = intrinsic_noise / len(her1) / 2
						
						# Fano factor: intrinsic noise * mean expression
						slice_fano_her1.append(intrinsic_noise * her1_mean)
						slice_fano_her7.append(intrinsic_noise * her7_mean)					
						
	# Bin Fano factor
	[binned_fixednumslice_mean_her1, binned_fixednumslice_fano_her1] = bin_fixednumslice(slice_mean_her1, slice_fano_her1)
	[binned_fixednumslice_mean_her7, binned_fixednumslice_fano_her7] = bin_fixednumslice(slice_mean_her7, slice_fano_her7)
	
	# Plot
	fig = plt.figure(figsize=(6,5),dpi=300)
	ax = fig.add_subplot(111)
	[her1,xmax_her1,ymax_her1] = plot(ax,binned_fixednumslice_mean_her1, binned_fixednumslice_fano_her1,'b') # circle
	[her7,xmax_her7,ymax_her7] = plot(ax,binned_fixednumslice_mean_her7, binned_fixednumslice_fano_her7,'r') # square
	
	xmax = max(xmax_her1,xmax_her7)
	ymax = max(ymax_her1,ymax_her7)	
	ax.xaxis.set_ticks(numpy.arange(0, xmax*1.1, determineTickInterval(xmax,5)))
	ax.set_xlim([0,xmax*1.05])
	ax.yaxis.set_ticks(numpy.arange(0, ymax*1.1, determineTickInterval(ymax,5)))
	ax.set_ylim([0,ymax*1.05])
		
	ax.legend([her1,her7],[r'\textit{her1}',r'\textit{her7}'], loc=2, scatterpoints=1, fontsize=12)
	ax.set_xlabel(r"Mean mRNA levels")
	ax.set_ylabel(r"Fano factor (intrinsic noise $\times$ mean)")
	ax.tick_params(direction='in')
	updateTicklabels(ax)
	fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95,  wspace=None, hspace=.3)	
	fig.savefig(directory + "/fano_factor.png", format = "png", dpi=300)
	
	# Write		
	workbook = xlwt.Workbook(encoding="ascii")	
	write(workbook.add_sheet("Her1"), binned_fixednumslice_mean_her1, binned_fixednumslice_fano_her1)
	write(workbook.add_sheet("Her7"), binned_fixednumslice_mean_her7, binned_fixednumslice_fano_her7)
	workbook.save(directory + "/fano_factor.xls")
	
def usage():
	print ("plot_fano_factor.py: Invalid command-line arguments")
	print ("Format: python plot_fano_factor.py <number of embryos> <her1 background noise mean> <her7 background noise mean> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print ("Example: python plot_fano_factor.py 20 3.083 1.311 ../wildtypefulldataset/output/embryo1/slices.xls ../wildtypefulldataset/output/embryo2/slices.xls ... ../wildtypefulldataset/output/embryo20/slices.xls ../wildtypefulldataset/output")
	exit(1)

main()
