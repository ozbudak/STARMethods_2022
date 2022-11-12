"""
Plot mean expression vs. Fano factor in different genetic backgrounds
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
markers = ['o','s'] # circle and square
genes = [r'\textit{her1}',r'\textit{her7}']

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
	ylabels = [format(label, r',.0f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_fano_factor.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2] 
	    	
	if len(sys.argv)==3*num_geneticbackgrounds+3:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # fano_factor.xls files
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_fano_factor.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
	else:
		usage()
	
	# Set up figure
	f = plt.figure(figsize=(6,5), dpi=300)
	ax = f.add_subplot(111)	
	xmax = float('-inf')
	ymax = float('-inf')
			
	for i in range(num_geneticbackgrounds): # genetic background
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_fano_factor.py: Cannot open file "'+inputs[i]+'".')
			exit(1)			
		worksheet_names = workbook.sheet_names()
		
		for j in range(2): # gene (her1 or her7)
			current_worksheet = workbook.sheet_by_name(worksheet_names[j])	
			
			for k in range(5): # bins			
				row = list(current_worksheet.row(k+1))	
				
				# Plot mean expression vs. Fano factor				
				ax.scatter(row[2].value, row[4].value, s=30, marker=markers[j], edgecolors='none', c=colors[i])
				ax.errorbar(row[2].value, row[4].value, xerr=2*row[3].value, yerr=2*row[5].value, ls='none', c=colors[i], capsize=2)	
				
				xmax = max(xmax,row[2].value+2*row[3].value)
				ymax = max(ymax,row[4].value+2*row[5].value)	
					
	# Figure labels and axis ticks	
	ax.set_xlabel(r"Mean mRNA levels")
	ax.set_xlim([0,math.ceil(xmax/10)*10])
	ax.set_ylabel(r"Fano factor (intrinsic noise $\times$ mean)")
	ax.set_ylim([0,ymax*1.05])	
	updateTicklabels(ax)
	
	# Figure legend
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds): # genetic background
		for j in range(2): # gene (her1 or her7)
			legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker=markers[j], markerfacecolor=colors[i], markeredgecolor=colors[i]))
			legend_labels.append(r''+names[i]+' '+genes[j])
	plt.legend(legend_shapes, legend_labels, bbox_to_anchor=(0., 1., 1., 0), loc=3, fontsize=12, ncol=num_geneticbackgrounds, numpoints=1)
	
	# Save figure
	f.subplots_adjust(left=0.1, bottom=0.1, right=.95, top=.85,  wspace=None, hspace=0.3)
	f.savefig(output_directory + "/compare_fano_factor.png", format = "png", dpi=300)
		
def usage():
	print ("compare_fano_factor.py: Invalid command-line arguments")
	print ("Format: python compare_fano_factor.py <number of genetic backgrounds> <output directory> <fano_factor.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_fano_factor.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/fano_factor.xls ../deltacfulldataset/output/fano_factor.xls ../deltadfulldataset/output/fano_factor.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)

main()
