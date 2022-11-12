# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 19:19:19 2018

@author: My
"""

import os, sys
import pandas as pd
import numpy, math
import matplotlib.pyplot as plt
from scipy.stats import sem
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

# Determined ticks for log-log plots
xt = [math.log10(4),math.log10(6),math.log10(8),math.log10(10),math.log10(20),math.log10(40),math.log10(60),math.log10(80),math.log10(100),math.log10(200)]
yt = [math.log10(0.01),math.log10(0.05),math.log10(0.10),math.log10(0.20),math.log10(0.30),math.log10(0.60),math.log10(0.80)]

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

def updateLogTicklabels(ax): # notation needed for log-log plot
	xlabels = [r'$\textsf{%.0f}$' % 10**label for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(direction='in',axis='x', pad=10)
	ylabels = [r'$\textsf{%.2f}$' % 10**label for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)
	ax.tick_params(direction='in')

def figureLabels(axeslog): # set axis labels for log-log plot
	for i in range(3):
		axeslog[i].set_xlabel(r"Total \textit{her} mRNA",fontsize=25)
	
	axeslog[0].set_ylabel(r"Total noise",fontsize=25)
	axeslog[1].set_ylabel(r"Intrinsic noise",fontsize=25)
	axeslog[2].set_ylabel(r"Extrinsic noise",fontsize=25)

def figureAxisTicks(axeslog, xmax, ymax, xminflog, xmaxflog, yminflog, ymaxflog, xminlog, xmaxlog, yminlog, ymaxlog): # set axis ticks for all figures
	for i in range(3):
		axeslog[i].xaxis.set_ticks(xt)
		axeslog[i].set_xlim(xminlog-abs(xmaxlog-xminlog)*0.05, math.log10(200))
		axeslog[i].yaxis.set_ticks(yt)
		axeslog[i].set_ylim(yminlog[i]-abs(ymaxlog[i]-yminlog[i])*0.075, ymaxlog[i]+abs(ymaxlog[i]-yminlog[i])*0.075)
		if i == 1:
			axeslog[i].set_ylim(math.log10(0.01), ymaxlog[i]+abs(ymaxlog[i]-yminlog[i])*0.075)
		updateLogTicklabels(axeslog[i])

def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_CVsquareed_fillArea.py: Number of genetic backgrounds must be an integer.')
		exit(1)
	elif not shared.isInt(sys.argv[2]):
		print ('compare_CVsquareed_fillArea.py: Number of bins for noise plots must be an integer.')
		exit(1)
	
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine (currently 3, can be increased)
	num_bins = int(sys.argv[2])
	output_directory = sys.argv[3]
	if len(sys.argv)==3*num_geneticbackgrounds+4:
		inputs = sys.argv[4:num_geneticbackgrounds+4] # noise.xls files
		for f in inputs:
			if not os.path.isfile(f):
				print ("compare_noise.py: File '"+f+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+4:2*num_geneticbackgrounds+4] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+4:] # colors for plotting different genetic backgrounds
	else:
		usage()
    
	# Figures showing each noise level separately with her x-axis
	f1 = plt.figure(figsize=(3*1.75,4), dpi=300) # compare_total_noise_log.png
	f2 = plt.figure(figsize=(3*1.75,4), dpi=300) # compare_intrinsic_noise_log.png
	f3 = plt.figure(figsize=(3*1.75,4), dpi=300) # compare_extrinsic_noise_log.png
	figures = [f1, f2, f3]
	axeslog = [f1.add_subplot(111),f2.add_subplot(111),f3.add_subplot(111)]
	xminlog = float('inf')
	xmaxlog = float('-inf')
	yminlog = [float('inf')]*3
	ymaxlog = [float('-inf')]*3

	for i in range(num_geneticbackgrounds): # genetic background
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_noise.py: Cannot open file "'+inputs[i]+'".')
			exit(1)
		worksheet_names = workbook.sheet_names()
			
		for j in range(4): # x-axis (her1, her7, her, harmonic mean)
			current_worksheet = workbook.sheet_by_name(worksheet_names[j])
		
			for k in range(num_bins): # bins
				# compare_noise.png
				row = list(current_worksheet.row(k+1))
				if row[1].value < 3:
					continue

		# compare_total_noise_log.png, compare_intrinsic_noise_log.png, compare_extrinsic_noise_log.png
		if j==2: # her x-axis
			for l in range(3): # noise
				axeslog[l].scatter(row[2].value, row[5+3*l].value,s = 22, edgecolors='none', c=colors[i])
				axeslog[l].plot(row[2].value, row[5+3*l].value)
				axeslog[l].errorbar(row[2].value, row[5+3*l].value,
									xerr=[[row[3].value],[row[4].value]],
									yerr=[[row[6+3*l].value],[row[7+3*l].value]], ls='none', c=colors[i],capsize=2, linewidth=0.5)
									yminlog[l] = min(yminlog[l],row[5+3*l].value-row[6+3*l].value)
										ymaxlog[l] = max(ymaxlog[l],row[5+3*l].value+row[7+3*l].value)
											xminlog = min(xminlog,row[2].value-row[3].value)
											xmaxlog = max(xmaxlog,row[2].value+row[4].value)

	# Figure axis labels, ticks, and legends
	figureLabels(axarr, axarrlog, axeslog)
	figureAxisTicks(axarr, axarrlog, axeslog, xmax, ymax, xminflog, xmaxflog, yminflog, ymaxflog, xminlog, xmaxlog, yminlog, ymaxlog)
	figureLegends(f, flog, num_geneticbackgrounds, names, colors)
		
	# Save figures
	f.subplots_adjust(left=0.05, bottom=0.05, right=.975, top=.95,  wspace=None, hspace=0.25)
	f.savefig(output_directory + "/compare_noise.png", format = "png", dpi=300)
		
	flog.subplots_adjust(left=0.05, bottom=0.05, right=.975, top=.95,  wspace=None, hspace=0.25)
	flog.savefig(output_directory + "/compare_noise_log.png", format = "png", dpi=300)
		
	#add legend for compare_*_noise.
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=colors[i], markeredgecolor=colors[i]))
		legend_labels.append(r''+names[i])
	
	for figure in figures:
		figure.legend(legend_shapes, legend_labels, bbox_to_anchor=(0., 1., 1., 0), loc=3, fontsize=12, ncol=num_geneticbackgrounds, numpoints=1)
		figure.subplots_adjust(left=0.14, bottom=0.15, right=.95, top=.85,  wspace=None, hspace=None)
	f1.savefig(output_directory + "/compare_total_noise_log.png", format = "png", dpi=300)
	f2.savefig(output_directory + "/compare_intrinsic_noise_log.png", format = "png", dpi=300)
	f3.savefig(output_directory + "/compare_extrinsic_noise_log.png", format = "png", dpi=300)
        
    
def usage():
	print ("plot_CVsquared_fillArea.py: Invalid command-line arguments")
	print ("Format: python plot_CVsquared_filled.py <output_directory> <list of colors to be used>")
	print ("Example: python compare_CVsquared.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/CVsquared.xls ../deltacfulldataset/output/CVsquared.xls ../deltadfulldataset/output/CVsquared.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)
    
main()
