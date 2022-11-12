"""
Plot mean expression vs. noise in different genetic backgrounds
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
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import sem
import xlrd
from xlrd import XLRDError
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

# Determined ticks for log-log plots
xt = [math.log10(4),math.log10(6),math.log10(8),math.log10(10),math.log10(20),math.log10(40),math.log10(60),math.log10(80),math.log10(100),math.log10(200)]
yt = [math.log10(0.01),math.log10(0.05),math.log10(0.10),math.log10(0.20),math.log10(0.30),math.log10(0.60),math.log10(0.80)]

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2,5,10,20,25,50]
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
	ax.tick_params(direction='in')
	
def updateLogTicklabels(ax): # notation needed for log-log plot
	xlabels = [r'$\textsf{%.0f}$' % 10**label for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(direction='in',axis='x', pad=10)
	ylabels = [r'$\textsf{%.2f}$' % 10**label for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)
	ax.tick_params(direction='in')
	
def figureLabels(axarr, axarrlog, axeslog): # set axis labels for all figures
	for i in range(3):
		axarr[i][0].set_xlabel(r"Total \textit{her1} RNA",fontsize=20)
		axarr[i][1].set_xlabel(r"Total \textit{her7} RNA",fontsize=20)
		axarr[i][2].set_xlabel(r"Total \textit{her} RNA",fontsize=20)
		axarr[i][3].set_xlabel(r"Total harmonic mean of \textit{her1} and \textit{her7}",fontsize=20)
		axarr[0][0].set_ylabel(r"Total noise",fontsize=20)
		axarr[1][0].set_ylabel(r"Intrinsic noise",fontsize=20)
		axarr[2][0].set_ylabel(r"Extrinsic noise",fontsize=20)		
		
		axarrlog[i][0].set_xlabel(r"Total \textit{her1} RNA",fontsize=20)
		axarrlog[i][1].set_xlabel(r"Total \textit{her7} RNA",fontsize=20)
		axarrlog[i][2].set_xlabel(r"Total \textit{her} RNA",fontsize=20)
		axarrlog[i][3].set_xlabel(r"Total harmonic mean of \textit{her1} and \textit{her7}",fontsize=20)
		axarrlog[0][0].set_ylabel(r"Total noise",fontsize=25)
		axarrlog[1][0].set_ylabel(r"Intrinsic noise",fontsize=25)
		axarrlog[2][0].set_ylabel(r"Extrinsic noise",fontsize=25)	
			
		axeslog[i].set_xlabel(r"Total \textit{her} mRNA",fontsize=25)
		
	axeslog[0].set_ylabel(r"Total noise",fontsize=25)
	axeslog[1].set_ylabel(r"Intrinsic noise",fontsize=25)
	axeslog[2].set_ylabel(r"Extrinsic noise",fontsize=25)

def figureAxisTicks(axarr, axarrlog, axeslog, xmax, ymax, xminflog, xmaxflog, yminflog, ymaxflog, xminlog, xmaxlog, yminlog, ymaxlog): # set axis ticks for all figures
	for i in range(3):			
		for j in range(4):
			axarr[i][j].xaxis.set_ticks(numpy.arange(0,xmax[j]*1.1,determineTickInterval(xmax[j],5)))
			axarr[i][j].set_xlim(0,xmax[j]*1.05)			
			axarr[i][j].yaxis.set_ticks(numpy.arange(0,ymax[i*4+j]*1.1,determineTickInterval(ymax[i*4+j],5)))
			axarr[i][j].set_ylim(0,ymax[i*4+j]*1.05)
			updateTicklabels(axarr[i][j])
			axarrlog[i][j].xaxis.set_ticks(xt)
			axarrlog[i][j].set_xlim(xminflog[j]-abs(xmaxflog[j]-xminflog[j])*0.05, xmaxflog[j]+abs(xmaxflog[j]-xminflog[j])*0.05)
			axarrlog[i][j].yaxis.set_ticks(yt)
			axarrlog[i][j].set_ylim(yminflog[i*4+j]-abs(ymaxflog[i*4+j]-yminflog[i*4+j])*0.05, ymaxflog[i*4+j]+abs(ymaxflog[i*4+j]-yminflog[i*4+j])*0.05)
			updateLogTicklabels(axarrlog[i][j])	
		axeslog[i].xaxis.set_ticks(xt)
		axeslog[i].set_xlim(xminlog-abs(xmaxlog-xminlog)*0.05, math.log10(200))
		axeslog[i].yaxis.set_ticks(yt)
		axeslog[i].set_ylim(yminlog[i]-abs(ymaxlog[i]-yminlog[i])*0.075, ymaxlog[i]+abs(ymaxlog[i]-yminlog[i])*0.075)
		if i == 1:
			axeslog[i].set_ylim(math.log10(0.01), ymaxlog[i]+abs(ymaxlog[i]-yminlog[i])*0.075)
		updateLogTicklabels(axeslog[i])
        
def figureLegends(f, flog, num_geneticbackgrounds, names, colors): # set figure legend
	legend_shapes = []
	legend_labels = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=colors[i], markeredgecolor=colors[i]))
		legend_labels.append(r''+names[i])
	f.legend(legend_shapes, legend_labels, 'upper center', fontsize=12, ncol=num_geneticbackgrounds, numpoints=1)
	flog.legend(legend_shapes, legend_labels, 'upper center', fontsize=12, ncol=num_geneticbackgrounds, numpoints=1)
	
	
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_noise.py: Number of genetic backgrounds must be an integer.')
		exit(1)
	elif not shared.isInt(sys.argv[2]):
		print ('compare_noise.py: Number of bins for noise plots must be an integer.')
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
	
	# Figure showing all x-axes --compare_noise.png
	f, axarr = plt.subplots(3,4) 
	f.set_size_inches(18,12)
	xmax = [float('-inf')]*4
	ymax = [float('-inf')]*12 
	
	# Log-log version of f --compare_noise_log.png
	flog, axarrlog = plt.subplots(3,4)
	flog.set_size_inches(18,12)
	xminflog = [float('inf')]*4
	xmaxflog = [float('-inf')]*4
	yminflog = [float('inf')]*12
	ymaxflog = [float('-inf')]*12
	
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
			noise_arr_x = {0:[],1:[],2:[]} # 0-Total, 1-Intrinsic, 2-Extrinsic
			noise_arr_y = {0:[],1:[],2:[]}
			noise_arr_y_uperr = {0:[],1:[],2:[]}
			noise_arr_y_lowerr = {0:[],1:[],2:[]}
			for k in range(num_bins): # bins				
				# compare_noise.png	
				row = list(current_worksheet.row(k+1))
				if row[1].value < 3:
					continue
				
				for l in range(3): # noise (total, intrinsic, extrinsic)
					axarr[l][j].scatter(row[2].value, row[4+2*l].value, c=colors[i])
					axarr[l][j].errorbar(row[2].value, row[4+2*l].value, xerr=2*row[3].value, yerr=2*row[5+2*l].value, ls='none', c=colors[i],capsize=2, linewidth=0.5)
					xmax[j] = max(xmax[j],row[2].value+2*row[3].value)
					ymax[l*4+j] = max(ymax[l*4+j],row[4+2*l].value+2*row[5+2*l].value)
				
				# compare_noise_log.png
				row = list(current_worksheet.row(k+14))
				for l in range(3): # noise (logarithmic)
					axarrlog[l][j].scatter(row[2].value, row[5+3*l].value, s = 22, edgecolors='none', c=colors[i])
					axarrlog[l][j].errorbar(row[2].value, row[5+3*l].value,
						xerr=[[row[3].value],[row[4].value]], yerr=[[row[6+3*l].value],[row[7+3*l].value]], ls='none', c=colors[i],capsize=2, linewidth=0.5)
					yminflog[l*4+j] = min(yminflog[l*4+j],row[5+3*l].value-row[6+3*l].value)
					ymaxflog[l*4+j] = max(ymaxflog[l*4+j],row[5+3*l].value+row[7+3*l].value)								
				xminflog[j] = min(xminflog[j], row[2].value-row[3].value)
				xmaxflog[j] = max(xmaxflog[j], row[2].value+row[3].value)
				
				# compare_total_noise_log.png, compare_intrinsic_noise_log.png, compare_extrinsic_noise_log.png
				if j==2: # her x-axis
					for l in range(3): # noise
						noise_arr_x[l].append(row[2].value)
						noise_arr_y[l].append(row[5+3*l].value)
						noise_arr_y_lowerr[l].append(row[6+3*l].value)
						noise_arr_y_uperr[l].append(row[7+3*l].value)
						axeslog[l].scatter(row[2].value, row[5+3*l].value,s = 22, edgecolors='none', c=colors[i])
						axeslog[l].errorbar(row[2].value, row[5+3*l].value,
							xerr=[[row[3].value],[row[4].value]],
							yerr=[[row[6+3*l].value],[row[7+3*l].value]], ls='none', c=colors[i],capsize=2, linewidth=0.5)
						yminlog[l] = min(yminlog[l],row[5+3*l].value-row[6+3*l].value)
						ymaxlog[l] = max(ymaxlog[l],row[5+3*l].value+row[7+3*l].value)
					xminlog = min(xminlog,row[2].value-row[3].value)
					xmaxlog = max(xmaxlog,row[2].value+row[4].value)
			if j==2:
				#max_std = []
				#min_std = []
				# connect dots in compare_*_noise_log.png graph
				for l in range(3):
					#max_std[l] = noise_arr_y[l]._add_(noise_arr_y_uperr[l])
					#min_std[l] = noise_arr_y[l] - noise_arr_y_lowerr[l]
					axeslog[l].plot(noise_arr_x[l],noise_arr_y[l],c=colors[i])
						#axeslog[l].fill_between(noise_arr_x[l], max_std, min_std, where=max_std>=min_std, facecolor=color[i], alpha=0.25)
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
	print ("compare_noise.py: Invalid command-line arguments")
	print ("Format: python compare_noise.py <number of genetic backgrounds> <output_directory> <number of bins used in plot_noise.py> <noise.xls from each genetic backgrounds> <genetic background names> <list of colors to be used>")
	print ("Example: python compare_noise.py 3 5 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/noise.xls ../deltacfulldataset/output/noise.xls ../deltadfulldataset/output/noise.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)

main()
