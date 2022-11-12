"""
Plot spatial noise (cell position vs. CV^2) in different genetic backgrounds
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

input_worksheet_names = ["Her1", "Her7", "Her"]
spss_worksheet_names = ["spss_ANOVA_her1", "spss_ANOVA_her7", "spss_ANOVA_her"]

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
	ax.tick_params(axis='x', pad=5)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_CVsquared.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2]
	if len(sys.argv)==3*num_geneticbackgrounds+3:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # CVsquared.xls files
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_CVsquared.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
	else:
		usage()	
	
	# Set up figures
	f = plt.figure(figsize=(4.5,4.5),dpi=300) # figure using only her on the x-axis 
	ax_her = f.add_subplot(111)
	
	fig = plt.figure(figsize=(4.5,12)) # figure using three different x-axes (her1, her7, her)	
	axes = [fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)]
	xmax = float('-inf') # keep track of how long the x-axis should be
	ymax = [float('-inf')]*3 # keep track of how long the y-axis of each subplot should be
	
	# Set up output file
	write_workbook = xlwt.Workbook(encoding="ascii")	
	write_worksheets = [write_workbook.add_sheet("Her1_spss"), \
	write_workbook.add_sheet("Her7_spss"),write_workbook.add_sheet("Her_spss")]	
	plot_data_worksheets = [write_workbook.add_sheet("Her1_plot_data"), \
	write_workbook.add_sheet("Her7_plot_data"), write_workbook.add_sheet("Her_plot_data")]
	labels = ['Genetic background','Posterior CV^2','','Genetic background','Anterior CV^2']
	# write data for spss sheets
	for ws in write_worksheets:
		for i in range(len(labels)):	
			ws.write(0,i,labels[i])
	# write data for sheets that contain data used to create plots
	labels = ["Slide index"] + [(name + "_cv2") for name in names] + [(name + "_ste_cv2") for name in names]
	for ws in plot_data_worksheets:
		for i, label in enumerate(labels):
			ws.write(0, i, label)
			
	posterior_row_indices = [1, 1, 1] # row index for posterior columns in each worksheet [her1.her7/her]
	anterior_row_indices = [1, 1, 1] # row index for anterior columns in each worksheet [her1/her7/her]
	plot_data_row_indicies = [1,1,1] # row index to write plot data [her1/her7/her]
	
	current_file_len = 0
	for i in range(num_geneticbackgrounds): # genetic background
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_CVsquared.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		worksheets = workbook.sheets()
		for j, ws_name in enumerate(input_worksheet_names): # gene (her1, her7, her)
			worksheet = workbook.sheet_by_name(ws_name)
			file_len =  worksheet.nrows			
			xmax = max(xmax, file_len-1)
			
			for k in range(file_len - 1): # slice (cell position)
				#write the slice index
				if i == 0: # to avoid overwriting into excel , we only write slice index when we create 
							# plot for the first genetic background
					(plot_data_worksheets[j]).write(k + 1, 0, k)
					current_file_len = k
				elif k > current_file_len: # but if we realize we have more slices than in the first genetic background,
											# we have to write it out
					(plot_data_worksheets[j]).write(k + 1, 0, k)
				
				# get the row data
				row = list(worksheet.row(k+1))
				
				# Plot mean CV^2 and two standard errors
				axes[j].scatter(row[0].value, row[2].value, s = 22, edgecolors='none', c=colors[i])
				axes[j].errorbar(row[0].value, row[2].value, yerr=2*row[3].value, ls='none', c=colors[i], capsize=2)
				ymax[j] = max(ymax[j], row[2].value+2*row[3].value)		
				
				#write data
				(plot_data_worksheets[j]).write(k + 1, i + 1, row[2].value)
				(plot_data_worksheets[j]).write(k + 1, i + 1 + num_geneticbackgrounds, row[3].value)
	
				if j == 2: # "her" worksheet
					ax_her.scatter(row[0].value, row[2].value, s = 22, edgecolors='none', c=colors[i])
					ax_her.errorbar(row[0].value, row[2].value, yerr=2*row[3].value, ls='none', c=colors[i], capsize=2)
				
			# Write data for SPSS statistical analysis
			spss_worksheet = workbook.sheet_by_name(spss_worksheet_names[j])
			file_len = spss_worksheet.nrows		
			for k in range(file_len-1):
				row = list(spss_worksheet.row(k+1))				
				if row[0].value==1: # if posterior
					write_worksheets[j].write(posterior_row_indices[j], 0, i+1)
					write_worksheets[j].write(posterior_row_indices[j], 1, row[1].value)
					posterior_row_indices[j]+=1
				else: # if anterior
					write_worksheets[j].write(anterior_row_indices[j], 3, i+1)
					write_worksheets[j].write(anterior_row_indices[j], 4, row[1].value)
					anterior_row_indices[j]+=1		
									
	write_workbook.save(output_directory + "/compare_CVsquared.xls")		
		
	# Figure axis labels and ticks	
	axes[0].set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her1)}")
	axes[1].set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her7)}")	
	axes[2].set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her)}")
	for i in range(len(axes)):
		axes[i].set_xlabel(r"Cell position (posterior - anterior)")	
		axes[i].xaxis.set_ticks(numpy.arange(0,xmax+1,10))
		axes[i].set_xlim(-1,xmax)
		axes[i].yaxis.set_ticks(numpy.arange(0,ymax[i]+1,determineTickInterval(ymax[i],5)))
		axes[i].set_ylim(0,ymax[i]*1.1)
		updateTicklabels(axes[i])	
	
	ax_her.set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her)}")
	ax_her.set_xlabel(r"Cell position (posterior - anterior)")	
	ax_her.xaxis.set_ticks(numpy.arange(0,xmax+1,10))
	ax_her.set_xlim(-1,xmax)
	ax_her.yaxis.set_ticks(numpy.arange(0,ymax[2]+1,determineTickInterval(ymax[2],5)))
	ax_her.set_ylim(0,ymax[2]*1.1)
	updateTicklabels(ax_her)
		
	# Figure legend		
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=colors[i], markeredgecolor=colors[i]))
		legend_labels.append(r''+names[i])
	axes[0].legend(legend_shapes, legend_labels, numpoints=1, fontsize=12, loc=9, bbox_to_anchor=(0.5,1.15), ncol=num_geneticbackgrounds, borderaxespad=0)
	ax_her.legend(legend_shapes, legend_labels, numpoints=1, fontsize=12, loc=9, bbox_to_anchor=(0.5,1.1), ncol=num_geneticbackgrounds, borderaxespad=0)
	
	# Save figures		
	f.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=.9,  wspace=None, hspace=.3)	
	f.savefig(output_directory + "/compare_CVsquared_her.png", format = "png", dpi=300)	

	fig.subplots_adjust(left=0.15, bottom=0.05, right=0.95, top=.95,  wspace=None, hspace=.3)	
	fig.savefig(output_directory + "/compare_CVsquared.png", format = "png", dpi=300)	

def usage():
	print ("compare_CVsquared.py: Invalid command-line arguments")
	print ("Format: python compare_CVsquared.py <number of genetic backgrounds> <output_directory> <CVsquared.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_CVsquared.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/CVsquared.xls ../deltacfulldataset/output/CVsquared.xls ../deltadfulldataset/output/CVsquared.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)
		
main()
