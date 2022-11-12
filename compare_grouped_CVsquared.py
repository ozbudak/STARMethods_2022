"""
Plot spatial noise (cell position vs. CV^2) from low, medium, and high expression groups in different genetic backgrounds
Copyright (C) 2017 Ahmet Ay, Dong Mai, Soo Bin Kwon, Ha Vu

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
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

# in each genetic background's CVsquared.xls file, these are sheets' names that we will have ot take and 
# write data for SPSS check
raw_data_sheet_names = ["her_low", "her_medium", "her_high"]

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 0.1

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(direction='in',axis='x', pad=5)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def plot_and_write_mean_ste(num_geneticbackgrounds, output_directory, inputs, names, colors, write_ws):
	# Set up excel file to write out data of each genetic background
	labels = ["genetic background","slice index","low_exp cv2", "low_exp cv2 ste", \
	"med_exp cv2", "med_exp cv2 ste", "high_exp cv2", "high_exp cv2 ste"]
	for i, label in enumerate(labels):
		write_ws.write(0, i , label)
	row_index = 1
	# Set up figure	
	f, axarr = plt.subplots(1,3) # total 3 subplots	
	f.set_size_inches(12,5)
	xmax = float('-inf') # keep track of how long the x-axis should be
	ymax = [float('-inf')]*3 # keep track of how long the y-axis of each subplot should be
	
	for i in range(num_geneticbackgrounds): # genetic background
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_grouped_CVsquared.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		worksheets = workbook.sheets()
		
		worksheet = workbook.sheet_by_name("Her")
		file_len =  worksheet.nrows			
		xmax = max(xmax, file_len-1)
					
		for k in range(file_len-1): # slice (cell position)
			row = list(worksheet.row(k+1))
			write_ws.write(row_index, 0 , i + 1) # write genetic background
			write_ws.write(row_index, 1, k) # write slice index
			for l in range(3): # expression
				if row[10+6*l].value>0: 
					axarr[l].scatter(row[0].value, row[8+6*l].value, s = 22, edgecolors='none', c=colors[i])
					axarr[l].errorbar(row[0].value, row[8+6*l].value, yerr=2*row[9+6*l].value, ls='none', c=colors[i],capsize=2, elinewidth=0.5)					
					ymax[l] = max(ymax[l], row[8+6*l].value+2*row[9+6*l].value)	
					# write into excel file
					write_ws.write(row_index, 2+l*2, row[8+6*l].value)
					write_ws.write(row_index, 3+l*2, row[9+6*l].value)
			row_index += 1

	
	# Figure axis labels and ticks
	axarr[0].set_title('Low expression', fontsize=12)	
	axarr[1].set_title('Medium expression', fontsize=12)	
	axarr[2].set_title('High expression', fontsize=12)		
	for i in range(3): # column			
		axarr[i].set_ylabel(r"\textit{her} $\textsf{CV}^{\textsf{\small{2}}}$ (averaged across embryos)",fontsize=18)
		axarr[i].set_xlabel(r"Cell position (posterior - anterior)",fontsize=16)	
		axarr[i].xaxis.set_ticks(numpy.arange(0,xmax+1,10))
		axarr[i].set_xlim(-1,xmax)
		axarr[i].yaxis.set_ticks(numpy.arange(0,ymax[i]+1,determineTickInterval(ymax[i],6)))
		axarr[i].set_ylim(0,ymax[i]*1.1)
		updateTicklabels(axarr[i])
	
	# Figure legend	
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=colors[i], markeredgecolor=colors[i]))
		legend_labels.append(r''+names[i])
	f.legend(legend_shapes, legend_labels, 'upper center', fontsize=12, ncol=num_geneticbackgrounds, numpoints=1)
	
	f.subplots_adjust(left=0.075, bottom=0.1, right=.95, top=.85,  wspace=0.25, hspace=None)	
	f.savefig(output_directory + "/compare_grouped_CVsquared.png", format = "png", dpi=300)	

def write_spss (num_geneticbackgrounds, inputs, wb):
	# create 3 sheets for 3 levels of mean expression
	spss_worksheets = []
	labels = ["genetic background", "cv2"]
	for i, sheet_name in enumerate(raw_data_sheet_names):
		spss_worksheets.append(wb.add_sheet(sheet_name))
		#write labels
		for j, label in enumerate(labels):
			(spss_worksheets[i]).write(0, j, label)
			
	# for each level of mean her expression: low, medium, high
	for i, sheet_name in enumerate(raw_data_sheet_names):
		row_index = 1
		for j in range(num_geneticbackgrounds): # for each genetic background
			try:
				read_wb = xlrd.open_workbook(inputs[j], 'r')
			except XLRDError as e: 
				print ('compare_grouped_CVsquared.py: Cannot open file "'+inputs[i]+'".')
				exit(1)
			read_ws = read_wb.sheet_by_name(sheet_name) # get the sheet of low/medium/high level
			file_len =  read_ws.nrows			
			for k in range(1, file_len): # for each row in the thing
				row = list(read_ws.row(k)) # get the row data
				(spss_worksheets[i]).write(row_index, 0, j + 1) # write genetic background
				(spss_worksheets[i]).write(row_index, 1, row[1].value) # write cv2 value
				row_index += 1
				
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_grouped_CVsquared.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2]
	
	if len(sys.argv)==3*num_geneticbackgrounds+3:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # CVsquared.xls files
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_grouped_CVsquared.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
	else:
		usage()
		
	write_wb = xlwt.Workbook(encoding="ascii")	
	sum_ws = write_wb.add_sheet("Summary")
	# plot and write summary data
	plot_and_write_mean_ste(num_geneticbackgrounds, output_directory, inputs, names, colors, sum_ws)
	# write data for statistical tests
	write_spss(num_geneticbackgrounds, inputs, write_wb)
	write_wb.save(output_directory + "/compare_grouped_CVsquared_her.xls")

def usage():
	print ("compare_grouped_CVsquared.py: Invalid command-line arguments")
	print ("Format: python compare_grouped_CVsquared.py <number of genetic backgrounds> <output_directory> <CVsquared.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_grouped_CVsquared.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/CVsquared.xls ../deltacfulldataset/output/CVsquared.xls ../deltadfulldataset/output/CVsquared.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)
		
main()
