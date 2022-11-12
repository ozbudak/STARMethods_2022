"""
Plot a bar graph comparing average noise in different genetic backgrounds
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
import matplotlib.patches as mpatches
from xlrd import XLRDError
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering
width3 = 0.75 # bar width for plotting
width2 = 1 # bar width for plotting figures with 2 genetic backgrounds
width5 = 0.75 # bar width for plotting figures with 5 genetic backgrounds

in_index = 0 
ex_index = 1
tot_index = 2
in_input_index = 5 
ex_input_index = 6 
tot_input_index = 7

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 0.1
	
def updateYTicklabels(ax):
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def get_data_and_write_excel(standard_file, inputs, names, num_bg, ws_spss):
	raw_results = []
	nor_results = [] # 3D: [genetic background][what type of noise][slide]
	
	# write labels into the excel file
	labels = ["genetic background", "raw intrinsic", "raw extrinsic", "raw total", \
	"normalized intrinsic", "normalized extrinsic", "normalized total"]
	for i , label in enumerate(labels):
		ws_spss.write(0, i, label)

	# get average of the standard background first (DMSO or wildtype)
	try:
		std_wb = xlrd.open_workbook(standard_file,'r')
	except XLRDError as e:
		print ('compare_noise_bar.py: Cannot open file "'+standard_file+'".')
		exit(1)				
	std_ws = std_wb.sheet_by_name('Combined')	
	std_file_len =  std_ws.nrows
	standard_data = [[], [], []] # [in, ex, total][slices]
	for i in range(1, std_file_len): # slice
		row = list(std_ws.row(i))
		(standard_data[in_index]).append(float(row[in_input_index].value))
		(standard_data[ex_index]).append(float(row[ex_input_index].value))
		(standard_data[tot_index]).append(float(row[tot_input_index].value))
	std_means = []
	for i in range(3):
		std_means.append(numpy.mean(standard_data[i]))
		
	row_index = 1 # row index to write data
	# Now process for all the genetic backgrounds
	
	for i in range (num_bg): # each genetic background
		try:
			wb = xlrd.open_workbook(inputs[i])
		except XLRDError as e:
			print ("compare_noise_bar.py: Cannot open file "+ inputs[i] + ".")
			exit(1)
		ws = wb.sheet_by_name("Combined")
		file_len = ws.nrows
		raw_noises = [[], [], []] # [in, ex, tot][slice of all different embryos]
		nor_noises = [[], [], []] # [in, ex, tot][slice of all different embryos]
		for j in range(1, file_len):
			row = list(ws.row(j))
			intrinsic = float(row[in_input_index].value)
			extrinsic = float(row[ex_input_index].value)
			total = float(row[tot_input_index].value)
			(raw_noises[in_index]).append(intrinsic)
			(raw_noises[ex_index]).append(extrinsic)
			(raw_noises[tot_index]).append(total)
			inNor = intrinsic / std_means[in_index]
			exNor = extrinsic / std_means[ex_index]
			totNor = total / std_means[tot_index]
			(nor_noises[in_index]).append(inNor)
			(nor_noises[ex_index]).append(exNor)
			(nor_noises[tot_index]).append(totNor)
			#write data into excel to do statistical tests later
			line = [names[i], intrinsic, extrinsic, total, inNor, exNor, totNor]
			for k, line_data in enumerate(line):
				ws_spss.write(row_index, k , line_data)
			row_index += 1
		# save data into places for plotting later
		raw_results.append(raw_noises)
		nor_results.append(nor_noises)
	return raw_results, nor_results

def create_noise_bar_plot(data, names, colors, width, num_geneticbackgrounds, save_file_name, ylabel, ws):
	# Write the title line of the worksheet
	labels = ['Genetic background','Total noise mean','Total noise ste',\
	'Intrinsic noise mean','Intrinsic noise ste','Extrinsic noise mean', 'Extrinsic noise ste']
	for i in range(len(labels)):
		ws.write(0,i,labels[i])	
	# Create plot
	f = plt.figure(figsize=(num_geneticbackgrounds*(width + 1),4), dpi=300) # combined
	ax = f.add_subplot(111)
	y_max = float('-inf')
	for i in range(num_geneticbackgrounds):
		tot_mean = numpy.mean(data[i][tot_index])
		tot_yerr = 2*numpy.std(data[i][tot_index])/math.sqrt(len(data[i][tot_index]))
		in_mean = numpy.mean(data[i][in_index])
		in_yerr = 2*numpy.std(data[i][in_index])/math.sqrt(len(data[i][in_index]))
		ex_mean = numpy.mean(data[i][ex_index])
		ex_yerr = 2*numpy.std(data[i][ex_index])/math.sqrt(len(data[i][ex_index]))
		# write stuff
		line = [names[i], tot_mean, tot_yerr, in_mean, in_yerr, ex_mean, ex_yerr]
		for j, data_write in enumerate(line):
			ws.write(i + 1, j , data_write) 
			
		# Plot stuff
		ax.bar(i, tot_mean, width, yerr= tot_yerr, color=colors[i], ecolor='k',capsize=4)
		ax.bar(i+(num_geneticbackgrounds+1), in_mean, width, yerr= in_yerr, color=colors[i], ecolor='k',capsize=4) 
		ax.bar(i+2*(num_geneticbackgrounds+1), ex_mean, width, yerr= ex_yerr, color=colors[i], ecolor='k',capsize=4) 		
		# find the limit of the yaxis
		y_max = max([y_max, tot_mean + tot_yerr, in_mean + in_yerr, ex_mean + ex_yerr])
	# Axis labels and ticks
	m = math.floor(num_geneticbackgrounds/2)
	if num_geneticbackgrounds%2==1: # odd number of genetic backgrounds
		ax.set_xticks([m+width/2, m+(num_geneticbackgrounds+1)+width/2, m+2*(num_geneticbackgrounds+1)+width/2])
	else: # even number of genetic backgrounds
		n = m-(1-width)/2
		ax.set_xticks([n, n+(num_geneticbackgrounds+1), n+2*(num_geneticbackgrounds+1)])
	ax.set_xticklabels((r'\begin{center}Total \linebreak noise \end{center}', r'\begin{center}Intrinsic\linebreak noise\end{center}', r'\begin{center}Extrinsic\linebreak noise\end{center}'), fontsize = 16)
	ax.set_xlim(-width,3*(num_geneticbackgrounds+1)-0.5)
	ax.tick_params(direction='in',axis='x', which='both', bottom='off', pad=10)
	ax.set_yticks(numpy.arange(0, y_max+1, determineTickInterval(y_max,5)))
	ax.set_ylim(0,y_max*1.05)	
	ax.set_ylabel(ylabel, fontsize = 16)
	updateYTicklabels(ax)
			
	# figure legend
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(mpatches.Patch(color=colors[i]))
		legend_labels.append(r''+names[i])
	ax.legend(legend_shapes, legend_labels, loc=9, bbox_to_anchor=(0.5,1.15), ncol=num_geneticbackgrounds, numpoints=1, fontsize=14)
	# Save figure
	f.subplots_adjust(left=0.15, bottom=0.15, right=.95, top=.9,  wspace=None, hspace=None)
	f.savefig(save_file_name, format = "png", dpi=300)

def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_noise_bar.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2]
	if len(sys.argv)==3*num_geneticbackgrounds+4:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # noise.xls files     
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_noise_bar.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
		standard_inputs = sys.argv[-1] # the last input is the file name of the standard genetic background
	else:
		usage()	
	
	
	# Set up figure	
	if num_geneticbackgrounds == 2:
		width = width2
	elif num_geneticbackgrounds == 3:
		width = width3
	else:
		width = width5
	
	# Set up output Excel file
	write_workbook = xlwt.Workbook(encoding="ascii")
	raw_worksheet = write_workbook.add_sheet("Raw_data") # Raw data of slices' noises 
	nor_worksheet = write_workbook.add_sheet("Normalized_data")	# Normalized_data of slices' noises
	spss_worksheet = write_workbook.add_sheet("spss_data") # data used to run statistical tests to find p values and stuff
	
	# get the data
	raw_data , nor_data = get_data_and_write_excel(standard_inputs, inputs, names, num_geneticbackgrounds, spss_worksheet)
	# draw compare_raw_noise_bar.png and write into raw_worksheet
	create_noise_bar_plot(raw_data, names, colors, width, num_geneticbackgrounds, \
	output_directory + "/compare_raw_noise_bar.png", "Noise", raw_worksheet)
	# draw compare_nor_noise_bar.png and write into nor_worksheet
	create_noise_bar_plot(nor_data, names, colors, width, num_geneticbackgrounds, \
	output_directory + "/compare_nor_noise_bar.png", "Normalized Noise", nor_worksheet)
	
	write_workbook.save(output_directory + "/compare_noise_bar.xls")
	
def usage():
	print ("compare_noise_bar.py: Invalid command-line arguments")
	print ("Format: python compare_noise_bar.py <number of genetic backgrounds> <output_directory> <raw_noise.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_noise_bar.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/raw_noise.xls ../deltacfulldataset/output/raw_noise.xls ../deltadfulldataset/output/raw_noise.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)

main()
