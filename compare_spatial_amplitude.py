"""
Plot spatial amplitudes in different genetic backgrounds 
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
import matplotlib.patches as mpatches
import xlrd, xlwt
from xlrd import XLRDError
from matplotlib import rc # text style
rc('text', usetex=True) # activate latex text rendering

width3 = 0.60 # bar width for plotting
width2 = 0.75
width5 = 0.75
def updateTicklabels(ax):
	ylabels = [format(label, r'.0f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def plot(ax, ngb, mean_amp_data, mean_amp_std_data, names, colors, width): # plot bar graphs showing spatial amplitude in different genetic backgrounds
	for i in range(2): # gene (her1 and her7) 
		for j in range(ngb): # genetic background	
			ax.bar(i*(ngb+1)+j, mean_amp_data[i][j], width, yerr=2*mean_amp_std_data[i][j], color=colors[j], ecolor='k',capsize=4, linewidth=0.5)
			
	# Figure axis ticks and labels
	m = (ngb-1.0)/2 + width/2
	ax.set_xticks([m, m+(ngb+1)])
	ax.set_xticklabels((r'\textit{her1}', r'\textit{her7}'),fontsize = 16)
	ax.tick_params(axis='x', which='both', bottom='off', pad=10)	
	ax.set_xlim(0-width, 1*(ngb+1)+ngb)
	ax.set_ylabel('Spatial amplitude', fontsize = 16)
	ax.tick_params(axis='y',direction='in')
	updateTicklabels(ax)
	
	# Figure legend
	legend_shapes = []
	legend_labels = []
	for i in range(ngb):
		legend_shapes.append(mpatches.Patch(color=colors[i]))
		legend_labels.append(r''+names[i])
	ax.legend(legend_shapes, legend_labels, loc=9, bbox_to_anchor=(0.5,1.15), ncol=ngb, numpoints=1, fontsize=12)

def write(sheet, data, names): # write data to Excel
	sheet.write(0,0,"Genetic background")
	sheet.write(0,1,"Her1 amplitude")
	sheet.write(0,3,"Genetic background")
	sheet.write(0,4,"Her7 amplitude")
	row_index = 1
	for j in range(len(data[0])): # genetic background
		for k in range(len(data[0][j])):
			sheet.write(row_index, 0, j+1)
			sheet.write(row_index, 1, data[0][j][k]) # her1 amplitude			
			sheet.write(row_index, 3, j+1) 
			sheet.write(row_index, 4, data[1][j][k]) # her7 amplitdue			
			row_index += 1
	
def main():	
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_spatial_amplitude.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2] # the directory to put files about data that compare different genetic backgrounds
	
	if len(sys.argv)==3*num_geneticbackgrounds+3:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # noise.xls files
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_spatial_amplitude.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
	else:
		usage()	
	
	mean_amp_her1 = [] # Mean amplitude
	mean_amp_std_her1 = []
	mean_amp_her7 = []
	mean_amp_std_her7 = []	
	amp_her1 = [] # Amplitude measured every five cell positions
	amp_her7 = []
	
	# Set up output file
	write_workbook = xlwt.Workbook(encoding="ascii")
	write_worksheet = write_workbook.add_sheet("Sheet 1")	
	write_worksheet.write(0,0,'Her1 amplitude')
	write_worksheet.write(0,len(names)+1,'Her7 amplitude')
	for i in range(len(names)):
		write_worksheet.write(1,i,names[i])
		write_worksheet.write(1,len(names)+i+1,names[i])
	
	# Read data 
	for i in range(num_geneticbackgrounds): # genetic background	
		amp_her1.append([])
		amp_her7.append([])
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_spatial_amplitude.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
			 						
		worksheet = workbook.sheet_by_name('Her1') # her1 sheet	
		row = list(worksheet.row(0)) # first row
		for j in range(1,len(row)):
			if row[j].value=='Average': # record mean amplitude
				mean_amp_her1.append(worksheet.cell(1,j).value)
				mean_amp_std_her1.append(worksheet.cell(2,j).value)
			elif row[j].value!='': # record raw amplitude
				amp_her1[i].append(worksheet.cell(1,j).value)
				
		worksheet = workbook.sheet_by_name('Her7') # her7 sheet	
		row = list(worksheet.row(0)) # first row
		for j in range(1,len(row)):
			if row[j].value=='Average': # record mean amplitude
				mean_amp_her7.append(worksheet.cell(1,j).value)
				mean_amp_std_her7.append(worksheet.cell(2,j).value)	
			elif row[j].value!='': # record raw amplitude
				amp_her7[i].append(worksheet.cell(1,j).value)
		
	
	# Plot data	
	if num_geneticbackgrounds == 2:
		width = width2
	elif num_geneticbackgrounds == 3:
		width = width3
	else:
		width = width5
	f = plt.figure(figsize=(num_geneticbackgrounds*(width + 1) + 0.5,4), dpi=300)
	ax = f.add_subplot(111)	
	plot(ax, num_geneticbackgrounds, [mean_amp_her1,mean_amp_her7],[mean_amp_std_her1,mean_amp_std_her7], names, colors, width)
		
	# Save figure
	f.subplots_adjust(left=0.175, bottom=0.15, right=.975, top=.9,  wspace=None, hspace=None)
	f.savefig(output_directory + "/compare_spatial_amplitude.png", format = "png", dpi=300)
	
	# Write data for SPSS statistical analysis
	workbook = xlwt.Workbook(encoding="ascii")
	worksheet = workbook.add_sheet("Sheet 1")			
	write(worksheet, [amp_her1, amp_her7], names)			
	workbook.save(output_directory + "/compare_spatial_amplitude.xls")
	
def usage():
	print ("compare_spatial_amplitude.py: Invalid command-line arguments")
	print ("Format: python compare_spatial_amplitude.py <number of genetic backgrounds> <output directory> <spatial_amplitude.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_spatial_amplitude.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/spatial_amplitude.xls ../deltacfulldataset/output/spatial_amplitude.xls ../deltadfulldataset/output/spatial_amplitude.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)

main()
