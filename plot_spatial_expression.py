"""
Plot spatial expression of an embryo: cell position vs. RNA expression levels 
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
import xlrd, xlwt
from xlrd import XLRDError
from scipy.ndimage import gaussian_filter1d # Gaussian filtering package
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

colors = ['b','r','g'] # colors for her1, her7, and her
region_names = ['left', 'right']
left_index = 0
right_index = 1
h1_index = 0
h7_index = 1
her_index = 2
def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.1,0.5,1,2,5,10,20,30,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def smoothen(x): # apply Gaussian smoothing to a given vector
	window = 3 # window size
	sigma = 10 # sigma
	for idx in range(0, len(x) - window + 1):
		smooth_slice = gaussian_filter1d(x[idx:idx+window], sigma)
		for idxx in range(len(smooth_slice)):
			x[idx + idxx] = smooth_slice[idxx]	
	return x # return smoothened x vector	

def interpolate(x): # interpolate empty data points
	for i in range(1, len(x)-1): # looking for one empty data point
		if numpy.isnan(x[i]) and not numpy.isnan(x[i-1]) and not numpy.isnan(x[i+1]):
			x[i] = (x[i-1] + x[i+1])/2 # average of two outer neighbors
	for i in range(1, len(x)-2): # looking for two consecutive empty data points
		if numpy.isnan(x[i]) and numpy.isnan(x[i+1]) and not numpy.isnan(x[i-1]) and not numpy.isnan(x[i+2]):
			x[i] = (x[i-1] + x[i+2])/2 # average of two outer neighbors 
			x[i+1] = (x[i-1] + x[i+2])/2
	return x

def determineStartEnd(x): # determine where data starts/ends to be valid (mutant may start/end with slices with too few cells)
	i = 0
	while i < (len(x)) and numpy.isnan(x[i]):
		i+=1	
	j = len(x)-1
	while j > 0 and numpy.isnan(x[j]):
		j-=1	
	return i, j+1 # start, end

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()] # intergers
	ax.set_xticklabels(xlabels)
	updateYTicklabels(ax)
	
def updateYTicklabels(ax):
	ylabels = [format(label, r',.0f') for label in ax.get_yticks()] # intergers
	ax.set_yticklabels(ylabels)
	
def writeEmbryo(ws, region, embryo_her1, embryo_stderr_her1, embryo_her7, embryo_stderr_her7, embryo_her, embryo_stderr_her): # write data for an embryo
	labels = ["Cell position"] + ["Raw mean","Std error","Smooth mean"]*3			
	ws.write(0,0,"Left") if region==0 else ws.write(0,11,"Right")
		
	ws.write(0,1+region*11,"Her1")
	ws.write(0,4+region*11,"Her7")
	ws.write(0,7+region*11,"Her")
	for l in range(len(labels)):
		ws.write(1,l+region*11,labels[l])
	for l in range(len(embryo_her1)): # slice
		line = [l,embryo_her1[l],embryo_stderr_her1[l],"",embryo_her7[l],embryo_stderr_her7[l],"",embryo_her[l],embryo_stderr_her[l],""]
		if numpy.isnan(embryo_her1[l]): # if nan, we leave the data cell blank
			line[1:3] = ["",""]
		if numpy.isnan(embryo_her7[l]):
			line[4:6] = ["",""]
		if numpy.isnan(embryo_her[l]):
			line[7:9] = ["",""]
		for m in range(len(line)):
			ws.write(2+l,m+region*11,line[m])	

def plotRawEmbryo(ax, embryo_data, embryo_stderr, color): # plot raw expression for an embryo
	ax.scatter(range(len(embryo_data)), embryo_data, c=color, s = 22, edgecolors='none')
	for k in range(len(embryo_data)):
		ax.errorbar(k, embryo_data[k], yerr=2*embryo_stderr[k], ls='none', c=color,capsize=2)	
	start, end = determineStartEnd(embryo_data)	
	ax.set_xlim(-2+start,end+1)			
	ax.set_xlabel("Cell position (posterior - anterior)")	
	start, end = ax.get_ylim()					
	ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
	ax.set_ylim(-end*0.1,end)					
	updateTicklabels(ax)
					
def plotRawSmoothEmbryo(ax_raw, ax_smooth, embryo_data, embryo_stderr, color): # plot raw and smoothened expression data for an embryo
	# Plot raw expression data
	plotRawEmbryo(ax_raw, embryo_data, embryo_stderr, color)	
	
	# Smoothen expression data
	data_length = len(embryo_data)	
	embryo_data = interpolate(embryo_data) # fill in empty data points so that the data is continuous
	start, end = determineStartEnd(embryo_data)	
	smooth_embryo_data = smoothen(embryo_data[start:end])
	r = ['']*start+smooth_embryo_data+['']*(data_length-end) # needed for writing data to Excel
	
	# Plot smoothened expression data
	ax_smooth.scatter(range(start, end), smooth_embryo_data, c=color, s = 22, edgecolors='none')	
	
	# Labels, ticks
	ax_smooth.tick_params(direction='in')
	ax_smooth.set_xlim(-2+start,end+1)	
	ax_smooth.set_xlabel("Cell position (posterior - anterior)", family='sans-serif')		
	start, end = ax_raw.get_ylim()
	ax_smooth.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
	ax_smooth.set_ylim(-end*0.1,end)
	updateTicklabels(ax_smooth)
	
	return r # return smoothened data

def plotWriteEmbryo(axes_raw, axes_smooth, ws, region, embryo_her1, embryo_stderr_her1, embryo_her7, embryo_stderr_her7, embryo_her, embryo_stderr_her): # plot and write data from an embryo		
	# Write embryo data to Excel
	writeEmbryo(ws, region, embryo_her1, embryo_stderr_her1, embryo_her7, embryo_stderr_her7, embryo_her, embryo_stderr_her)
								
	# Plot raw and smoothened expression data		
	smooth_embryo_her1 = plotRawSmoothEmbryo(axes_raw[0], axes_smooth[0], embryo_her1, embryo_stderr_her1, 'b')
	smooth_embryo_her7 = plotRawSmoothEmbryo(axes_raw[1], axes_smooth[1], embryo_her7, embryo_stderr_her7, 'r')
	smooth_embryo_her = plotRawSmoothEmbryo(axes_raw[2], axes_smooth[2], embryo_her, embryo_stderr_her, 'g')
	axes_raw[0].set_ylabel(r"\textit{her1} mRNA")
	axes_raw[1].set_ylabel(r"\textit{her7} mRNA")
	axes_raw[2].set_ylabel(r"\textit{her} mRNA")		
	
	# Write smoothened data to Excel
	smooth_embryo_data = [smooth_embryo_her1, smooth_embryo_her7, smooth_embryo_her]
	for l in range(len(smooth_embryo_data)):
		for m in range(len(smooth_embryo_data[l])):
			ws.write(2+m,(l+1)*3+region*11,smooth_embryo_data[l][m])

def process_region(in_ws):
	'''
	in_ws: cells.xls
	'''
	###### 1. Process data #######
	mean = [[],[], []] #[h1,h7 or h][slice]
	stderr = [[],[], []] #[h1, h7 or h][slice]
	her1 = [] #[slice][cell]
	her7 = [] #[slice][cell]
	her = [] # [slice][cell]
	file_len = in_ws.nrows
	for k in range(1, file_len):
		row = list(in_ws.row(k)) # slice
		if isinstance(row[1].value,float): # enough cells to analyze
			num_cells = int(row[1].value) # number of cells in this slice
			slice_h1 = [] #[cell]
			slice_h7 = []
			slice_her = []
			
			for l in range(num_cells):
				h1 = row[8+2*l].value # cell's her1 expression level
				h7 = row[8+2*l+1].value # cell's her7 expression level
				
				# Set negative expression levels to zero (do not discard)
				slice_h1.append(h1 if h1 > 0 else 0) 
				slice_h7.append(h7 if h7 > 0 else 0)
				slice_her.append(h1+h7 if h1+h7 > 0 else 0)		
			#mean of slice
			mean[h1_index].append(numpy.mean(slice_h1))
			mean[h7_index].append(numpy.mean(slice_h7))
			mean[her_index].append(numpy.mean(slice_her))
			#ste of slice
			stderr[h1_index].append(numpy.std(slice_h1) / math.sqrt(num_cells))
			stderr[h7_index].append(numpy.std(slice_h7) / math.sqrt(num_cells))
			stderr[her_index].append(numpy.std(slice_her) / math.sqrt(num_cells))
			# all data of slice
			her1.append(slice_h1)
			her7.append(slice_h7)
			her.append(slice_her)
		else:
			#mean of slice
			mean[h1_index].append(numpy.nan)
			mean[h7_index].append(numpy.nan)
			mean[her_index].append(numpy.nan)
			#ste of slice
			stderr[h1_index].append(numpy.nan)
			stderr[h7_index].append(numpy.nan)
			stderr[her_index].append(numpy.nan)
			# all data of slice
			her1.append([])
			her7.append([])
			her.append([])
	return mean, stderr, her1, her7, her

def create_combined_figure(write_worksheet, directory, region_index, mean, stderr):
	'''
	create left_spatial_expression.xls and right_spatial_expression.xls
	'''
	figure_combined = plt.figure(figsize=(10,9), dpi=300)
	axes_raw = [figure_combined.add_subplot(321), figure_combined.add_subplot(323), figure_combined.add_subplot(325)]
	axes_smooth = [figure_combined.add_subplot(322), figure_combined.add_subplot(324), figure_combined.add_subplot(326)]							
	plotWriteEmbryo(axes_raw,axes_smooth,write_worksheet,region_index,mean[h1_index],stderr[h1_index],\
	mean[h7_index],stderr[h7_index],mean[her_index],stderr[her_index])
	figure_combined.savefig(directory + "/" + region_names[region_index] + "_spatial_expression.png", format="png", dpi = 300)
	
	
def create_plot_raw_embryo(directory, region_index, mean, stderr):
	figure_raw_her1 = plt.figure(figsize=(6.5,4), dpi=300)
	ax_raw_her1 = figure_raw_her1.add_subplot(111)
		
	figure_raw_her7 = plt.figure(figsize=(6.5,4), dpi=300)
	ax_raw_her7 = figure_raw_her7.add_subplot(111)
	
	plotRawEmbryo(ax_raw_her1, mean[h1_index], stderr[h1_index], colors[0])
	ax_raw_her1.set_ylabel(r"\textit{her1} mRNA")			
	figure_raw_her1.subplots_adjust(left=0.1, bottom=0.15, right=0.975, top=.95, wspace=None, hspace=None)		
	
	plotRawEmbryo(ax_raw_her7, mean[h7_index], stderr[h7_index], colors[1])
	ax_raw_her7.set_ylabel(r"\textit{her7} mRNA")			
	figure_raw_her7.subplots_adjust(left=0.1, bottom=0.15, right=0.975, top=.95, wspace=None, hspace=None)
	
	figure_raw_her1.savefig(directory + "/" + region_names[region_index] + "_raw_her1.png", format = "png", dpi=300)
	figure_raw_her7.savefig(directory + "/" + region_names[region_index] + "_raw_her7.png", format = "png", dpi=300)	

def plotScatterSlice(ax, her, color): # plot scatters of cells' her 1/7 mRNAs levels and the cells' position (slice)
	for i, slide in enumerate(her):
		ax.scatter([i] * len(slide), slide, c = color, s = 22, alpha = 0.2)
	ax.set_xlim(0-2, len(her) + 1)
	ax.set_xlabel("Cell position (posterior - anterior)")	
	start, end = ax.get_ylim()					
	ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
	ax.set_ylim(-end*0.1,end)					
	updateTicklabels(ax)

def create_scatter_slice(directory, region_index, her1, her7):
	figure_scat_her1 = plt.figure(figsize=(6.5,4), dpi=300)
	ax_h1 = figure_scat_her1.add_subplot(111)
	
	figure_scat_her7 = plt.figure(figsize=(6.5,4), dpi=300)
	ax_h7 = figure_scat_her7.add_subplot(111)
	
	plotScatterSlice(ax_h1, her1, colors[0])
	ax_h1.set_ylabel(r"\textit{her1} mRNA")
	figure_scat_her1.subplots_adjust(left=0.1, bottom=0.15, right=0.975, top=.95, wspace=None, hspace=None)		
	
	plotScatterSlice(ax_h7, her7, colors[1])
	ax_h7.set_ylabel(r"\textit{her7} mRNA")
	figure_scat_her7.subplots_adjust(left=0.1, bottom=0.15, right=0.975, top=.95, wspace=None, hspace=None)		
	
	figure_scat_her1.savefig(directory + "/" + region_names[region_index] + "_scat_her1.png", format = "png", dpi=300)
	figure_scat_her7.savefig(directory + "/" + region_names[region_index] + "_scat_her7.png", format = "png", dpi=300)	

def main():
	args = sys.argv[1:]
	num_args = len(args)
	req_args = [False]*2
	if num_args == 4:
		for arg in range(0, num_args - 1, 2):
			option = args[arg]
			value = args[arg + 1]			
			if option == '-i' or option == '--input-file': # slices.xls
				filename = value
				if not os.path.isfile(filename):
					print("plot_spatial_expression.py: File "+filename+" does not exist.")
					exit(1)
				req_args[0] = True
			elif option == '-d' or option == '--output-directory':
				directory = value
				req_args[1] = True							
			elif option == '-h' or option == '--help':
				usage()
			else:
				usage()
		for arg in req_args:
			if not arg:
				usage()
	else:
		usage()
	
	shared.ensureDir(directory)
	
	try:
		workbook = xlrd.open_workbook(filename,'r')
	except XLRDError as e:
		print("plot_spatial_expression.py: Cannot open"+filename+".")
		exit(1)	
	worksheets = workbook.sheets()	
	
	# Set up output Excel file
	write_workbook = xlwt.Workbook(encoding="ascii") # output Excel file
	write_worksheet = write_workbook.add_sheet("Sheet 1", cell_overwrite_ok=True)
	# For each region
	for j in range(len(worksheets)): # region left of right 
		worksheet = worksheets[j]
		# get input
		mean, stderr, her1, her7, her = process_region(worksheet)
		# create left or right _spatial_expression.png
		create_combined_figure(write_worksheet, directory, j, mean, stderr)
		# create left or right_raw_her1.png and ..._raw_her7.png
		create_plot_raw_embryo(directory, j, mean, stderr)
		# create left or right _scat_her1.png and ..._scat_her7.png, 
		# we comment this out because the figures do not add much new information to the study
		#create_scatter_slice(directory, j, her1, her7)
		
	write_workbook.save(directory+"/spatial_expression.xls")
	
def usage():
	print("plot_spatial_expression.py: Invalid command-line arguments")
	print("Format: python plot_spatial_expression.py -i <input Excel file> -d <output directory")
	print("Example: python plot_spatial_expression -i ../wildtypefulldataset/output/embryo1/slices.xls -d ../wildtypefulldataset/output/embryo1")
	exit(1)

main()
