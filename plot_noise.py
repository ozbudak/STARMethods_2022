"""
Plot mean expression vs. noise 
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
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.legend_handler as handler
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import xlrd, xlwt
from xlrd import XLRDError
from matplotlib import rc # text style 
from decimal import Decimal
rc('text', usetex=True) # activate latex text rendering
colors = ['#2A6EFF','m','g', 'none']
bin_range = [45,80,100,125,300] # Fix bin range defined
#bin_range = [20,30,55,70,300] # Fix bin range for nrarp 
#bin_range = [20,30,300]#Fix bin range for nrarp wt vs. dletaC
#bin_range = [20,40,60,150] # Fix bin range for wt her1v her7 single gene noise 
#bin_range = [20,40,60,80,150,300] # Fix bin range for deltaC, her single gene noise
#background_bin_tb = [0.071901823, 0.059789701, 0.054201522, 0.059651073, 0.057268856] #total measurement error for no volume correction
#background_bin_tb = [0.031906776, 0.030027575, 0.030061986, 0.033075336, 0.040620712] #total measurement error for volume correction
#background_bin_im = [0.03, 0.03, 0.03, 0.03, 0.03] # intrinsic measurement error
#background_bin_eg_novolume = [0.041901823,0.029789701, 0.024201522, 0.029651073, 0.027268856] # extrinsic measurement error for no volume correction
#background_bin_eg = [0.001906776,0.0000275752, 0.0000619856, 0.003075336, 0.010620712] # extrinsic measurement error for volume correction

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [1,2,5,10,20,30,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def binData(slice_mean, slice_intrinsic, slice_extrinsic, slice_total, num_bins): # group noise levels based equal slices bin
	bin_size = (max(slice_mean)-min(slice_mean))/num_bins		
	binned_x = []
	binned_intrinsic = []
	binned_extrinsic = []
	binned_total = []	
	for i in range(num_bins):
		binned_x.append([])
		binned_intrinsic.append([])
		binned_extrinsic.append([])
		binned_total.append([])		
	for i in range(len(slice_mean)):
		for j in range(num_bins):
			if slice_mean[i]<=min(slice_mean)+bin_size*(j+1):
				binned_x[j].append(slice_mean[i])
				binned_intrinsic[j].append(slice_intrinsic[i])
				binned_extrinsic[j].append(slice_extrinsic[i])
				binned_total[j].append(slice_total[i])				
				break						
	return 	binned_x, binned_intrinsic, binned_extrinsic, binned_total

def binData_fix(slice_mean, slice_intrinsic, slice_extrinsic, slice_total): # group noise levels based on Fixed expression value
	num_bins = len(bin_range)
	binned_x = []
	binned_intrinsic = []
	binned_extrinsic = []
	binned_total = []	
	for i in range(num_bins):
		binned_x.append([])
		binned_intrinsic.append([])
		binned_extrinsic.append([])
		binned_total.append([])		
	for i in range(len(slice_mean)):
		for j in range(num_bins):
			if slice_mean[i]<=bin_range[j]:
				binned_x[j].append(slice_mean[i])
				binned_intrinsic[j].append(slice_intrinsic[i])
				binned_extrinsic[j].append(slice_extrinsic[i])
				binned_total[j].append(slice_total[i])				
				break						
	return	binned_x, binned_intrinsic, binned_extrinsic, binned_total

def binData_fix_single_her1_gene(slice_mean, slice_total): # group noise levels based on Fixed expression value
	num_bins = len(bin_range)
	binned_x = []
	binned_total = []	
	for i in range(num_bins):
		binned_x.append([])
		binned_total.append([])		
	for i in range(len(slice_mean)):
		for j in range(num_bins):
			if slice_mean[i]<=bin_range[j]:
				binned_x[j].append(slice_mean[i])
				binned_total[j].append(slice_total[i])				
				break						
	return	binned_x, binned_total

def binData_single_gene_fix(slice_mean, slice_total): # Oriana group single gene noise levels based on Fixed single gene expression value 
	num_bins = len(bin_range)
	binned_x = []
	binned_total = []	
	for i in range(num_bins):
		binned_x.append([])
		binned_total.append([])	
	for i in range(len(slice_mean)):
		for j in range(num_bins):
			if slice_mean[i]<=bin_range[j]:
				binned_x[j].append(slice_mean[i])
				binned_total[j].append(slice_total[i])				
				break						
	return	binned_x,binned_total

def writeSPSS5_bin(ws, binned_x, binned_intrinsic, binned_extrinsic, binned_total, binsize): # write data for SPSS statistical
	for n in range(binsize):
		labels = ['#Bin'+str(n+1)+'_SliceMean','Total noise level','Intrinsic noise level','Extrinsic noise level','']
		for i in range(len(labels)):
			ws.write(0,i+n*5,labels[i])

		for j in range(len(binned_x[n])):
			ws.write(j+1,0+n*5,binned_x[n][j])
			ws.write(j+1,1+n*5,binned_total[n][j])
			ws.write(j+1,2+n*5,binned_intrinsic[n][j])
			ws.write(j+1,3+n*5,binned_extrinsic[n][j])		

def writeSPSS5_single_bin(ws, binned_x, binned_total, binsize): # Oriana write single gene noise data for SPSS statistical
	for n in range(binsize):
		labels = ['#Bin'+str(n+1)+'_SliceMean','Total noise level','']
		for i in range(len(labels)):
			ws.write(0,i+n*5,labels[i])

		for j in range(len(binned_x[n])):
			ws.write(j+1,0+n*5,binned_x[n][j])
			ws.write(j+1,1+n*5,binned_total[n][j])

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(direction='in',axis='x', pad=10,labelsize = 16)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def updateLogTicklabels(ax): # notation needed for log-log plot
	xlabels = [r'$\textbf{%.0f}$' % 10**label for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=10)
	ylabels = [r'$\textbf{%.2f}$' % 10**label for label in ax.get_yticks()] 
	ax.set_yticklabels(ylabels)

def findLogXlim(xmin, xmax): # find appropriate limits for logarithmic x-axis given its minimum and maximum
	xlim = [math.log10(4), math.log10(400)] # largest possible range
	# Adjust x-axis lower limit (left end) based on data
	if xmin>math.log10(20)+0.02:
		xlim[0] = math.log10(20)	# 20
	elif xmin>math.log10(15)+0.02:
		xlim[0] = math.log10(15)	# 15
	elif xmin>1.02:
		xlim[0] = 1			# 10
	elif xmin>math.log10(8)+0.02:
		xlim[0] = math.log10(8)		# 8
	elif xmin>math.log10(6)+0.02:
		xlim[0] = math.log10(6)		# 6
	elif xmin>math.log10(4)+0.02:
		xlim[0] = math.log10(4)		# 4
	elif xmin>math.log10(2)+0.02:
		xlim[0] = math.log10(2)		# 2		
	# Adjust x-axis upper limit (right end) based on data
	if xmax<math.log10(20)-0.02:
		xlim[1] = math.log10(20)	# 20 
	elif xmax<math.log10(40)-0.02:
		xlim[1] = math.log10(40)	# 40
	elif xmax<math.log10(60)-0.02:
		xlim[1] = math.log10(60)	# 60
	elif xmax<math.log10(80)-0.02:
		xlim[1] = math.log10(80)	# 80
	elif xmax<2-0.02:
		xlim[1] = 2			# 100
	elif xmax<math.log10(150)-0.02:
		xlim[1] = math.log10(150)	# 150
	elif xmax<math.log10(200)-0.02:
		xlim[1] = math.log10(200)	# 200
	return xlim
	
def plotWrite(sheet, binned_x, binned_intrinsic, binned_extrinsic, binned_total, xlabel, directory): # plot and write mean expression vs. noise	
	fig = plt.figure(figsize = (12,7),dpi=300)
	ax = fig.add_subplot(111)
	xmin = float('inf')
	xmax = float('-inf')
	ymax = float('-inf')
	
	mean_level = []
	total_noise = []
	intrinsic_noise = []
	extrinsic_noise = []
	x_error = []
	y_tot_error = []
	y_in_error = []
	y_ex_error = []
	
	# Calculate ys
	for i in range(len(binned_x)):
		xmin = min(xmin,numpy.mean(binned_x[i]))
		xmax = max(xmax,numpy.mean(binned_x[i])+2*numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])))
		ymax = max(ymax,numpy.mean(binned_total[i])+2*numpy.std(binned_total[i])/math.sqrt(len(binned_total[i])))
		mean_level.append(numpy.mean(binned_x[i]))
		total_noise.append(numpy.mean(binned_total[i]))
		intrinsic_noise.append(numpy.mean(binned_intrinsic[i]))
		extrinsic_noise.append(numpy.mean(binned_extrinsic[i]))
		x_error.append(2*numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])))
		y_tot_error.append(2*numpy.std(binned_total[i])/math.sqrt(len(binned_total[i])))
		y_in_error.append(2*numpy.std(binned_intrinsic[i])/math.sqrt(len(binned_intrinsic[i])))
		y_ex_error.append(2*numpy.std(binned_extrinsic[i])/math.sqrt(len(binned_extrinsic[i])))
	# Plot
	
	tot = ax.scatter(mean_level, total_noise, s = 22, edgecolors='none', c=colors[0], label = "Total")
	intr = ax.scatter(mean_level, intrinsic_noise, s = 22, edgecolors='none', c=colors[1], label = "Intrinsic")
	extr = ax.scatter(mean_level, extrinsic_noise, s = 22, edgecolors='none', c=colors[2], label = "Extrinsic")
	ax.legend(loc = 0, scatterpoints=1, fontsize = 14)
	'''if (len(binned_x) == len(background_bin_tb)):
		background_tb = ax.plot(mean_level, background_bin_tb, '--', c = colors[0], label="Baseline Noise")
		background_eg = ax.plot(mean_level, background_bin_im, '--', c = colors[1], label="Baseline Noise")
		background_im = ax.plot(mean_level, background_bin_eg, '--', c = colors[2], label="Baseline Noise")'''
	
	ax.errorbar(mean_level, total_noise, xerr = x_error, yerr = y_tot_error, ls='none', c=colors[0],capsize=2, elinewidth=0.5)
	ax.errorbar(mean_level, intrinsic_noise, xerr = x_error, yerr = y_in_error, ls='none', c=colors[1],capsize=2, elinewidth=0.5)
	ax.errorbar(mean_level, extrinsic_noise, xerr = x_error, yerr = y_ex_error, ls='none', c=colors[2],capsize=2, elinewidth=0.5)
	'''if len(binned_x) == len(background_bin_tb):
		ax.scatter(mean_level, background_bin_tb, s = 22, edgecolors='none', c=colors[0])
		ax.scatter(mean_level, background_bin_im, s = 22, edgecolors='none', c=colors[1])
		ax.scatter(mean_level, background_bin_eg, s = 22, edgecolors='none', c=colors[2])'''

	ax.tick_params(direction='in',labelsize = 14)
	ax.xaxis.set_ticks(numpy.arange(0, xmax+10, determineTickInterval(xmax,10)))
	if xmax>=1: # need more space on the x-axis
		ax.set_xlim(0,xmax+xmin)
	else:
		ax.set_xlim(0,1)
	ax.yaxis.set_ticks(numpy.arange(0, ymax+.1, .2), determineTickInterval(ymax+ 0.1, 10))
	ax.set_ylim(0,ymax+0.05)
	ax.set_ylabel("Noise (coefficient of variation squared)", fontsize = 17)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in grouped slices" % {'label':xlabel}, fontsize = 10)

	# save figure
	fig.subplots_adjust(left=0.125, bottom=0.15, right=.95, top=.95, wspace=0.2, hspace=0.3)
	fig.savefig(directory + "/noise_" + xlabel + ".png", format = "png", dpi=300)
	# Write data
	labels = ['Bin','# slices','Normalized mean','Std error','Total noise','Std error','Intrinsic noise','Std error','Extrinsic noise','Std error']
	for i in range(len(labels)):
		sheet.write(0,i,labels[i])
	for i in range(len(binned_x)):
		line = [i+1, len(binned_x[i]), numpy.mean(binned_x[i]), numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])),
			numpy.mean(binned_total[i]), numpy.std(binned_total[i])/math.sqrt(len(binned_total[i])),
			numpy.mean(binned_intrinsic[i]), numpy.std(binned_intrinsic[i])/math.sqrt(len(binned_intrinsic[i])),
			numpy.mean(binned_extrinsic[i]), numpy.std(binned_extrinsic[i])/math.sqrt(len(binned_extrinsic[i]))]
		for j in range(len(line)):
			sheet.write(i+1,j,line[j])	
				

def plotLog(ax, binned_x, binned_intrinsic, binned_extrinsic, binned_total): # plot mean expression vs. noise in logarithmic scale (no writing)	
	xmin = float('inf') 
	xmax = float('-inf')
	ymin = float('inf') 
	ymax = float('-inf')
		
	# Plot data
	for i in range(len(binned_x)):	
		x = math.log10(numpy.mean(binned_x[i]))
		total = math.log10(numpy.mean(binned_total[i]))
		intrinsic = math.log10(numpy.mean(binned_intrinsic[i]))
		extrinsic = math.log10(numpy.mean(binned_extrinsic[i]))
		
		# Standard errors
		xe = numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])) 
		total_ye = numpy.std(binned_total[i])/math.sqrt(len(binned_total[i]))
		intrinsic_ye = numpy.std(binned_intrinsic[i])/math.sqrt(len(binned_intrinsic[i])) 
		extrinsic_ye = numpy.std(binned_extrinsic[i])/math.sqrt(len(binned_extrinsic[i]))
		
		# Applying log to 2 standard errors
		# Transforming error bars in log scale: http://labs.physics.dur.ac.uk/skills/skills/logscales.php
		xe_left = abs(math.log10(abs(numpy.mean(binned_x[i])-2*xe))-x) 
		xe_right = abs(math.log10(abs(numpy.mean(binned_x[i])+2*xe))-x)
		total_ye_lower = abs(math.log10(abs(numpy.mean(binned_total[i])-2*total_ye))-total)
		total_ye_upper = abs(math.log10(abs(numpy.mean(binned_total[i])+2*total_ye))-total)
		intrinsic_ye_lower = abs(math.log10(abs(numpy.mean(binned_intrinsic[i])-2*intrinsic_ye))-intrinsic)
		intrinsic_ye_upper = abs(math.log10(abs(numpy.mean(binned_intrinsic[i])+2*intrinsic_ye))-intrinsic)
		extrinsic_ye_lower = abs(math.log10(abs(numpy.mean(binned_extrinsic[i])-2*extrinsic_ye))-extrinsic)
		extrinsic_ye_upper = abs(math.log10(abs(numpy.mean(binned_extrinsic[i])+2*extrinsic_ye))-extrinsic)
		
		# Find minimum and maximum of data points including error bars 				
		xmin = min(xmin, x-xe_left)
		xmax = max(xmax, x+xe_right)
		ymin = min(ymin, intrinsic-intrinsic_ye_lower, extrinsic-extrinsic_ye_lower)
		ymax = max(ymax, total+total_ye_upper)
		
		# Total
		ax.scatter(x, total, s = 22, edgecolors='none', c=colors[0])
		ax.errorbar(x, total, xerr=[[xe_left],[xe_right]], yerr=[[total_ye_lower],[total_ye_upper]], ls='none', c=colors[0],capsize=2, elinewidth=0.5)	
		# Intrinsic
		ax.scatter(x, intrinsic, s = 22, edgecolors='none', c=colors[1])
		ax.errorbar(x, intrinsic, xerr=[[xe_left],[xe_right]], yerr=[[intrinsic_ye_lower],[intrinsic_ye_upper]], ls='none', c=colors[1],capsize=2, elinewidth=0.5)
		# Extrinsic
		ax.scatter(x, extrinsic, s = 22, edgecolors='none', c=colors[2])
		ax.errorbar(x, extrinsic, xerr=[[xe_left],[xe_right]], yerr=[[extrinsic_ye_lower],[extrinsic_ye_upper]], ls='none', c=colors[2],capsize=2, elinewidth=0.5)
	
	xt = [1,2,4,6,8,10,20,40,60,80,100,200,400]
	for i in range(len(xt)):
		xt[i] = math.log10(xt[i])	
	ax.xaxis.set_ticks(xt)
	ax.set_xlim(findLogXlim(xmin,xmax)) # use minimum and maximum to find appropriate x-axis limits
		
	yt = [0.01, 0.001,0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.00]
	for i in range(len(yt)):
		yt[i] = math.log10(yt[i])
	ax.yaxis.set_ticks(yt)
	ax.set_ylim([-2, ymax+0.1])	
	ax.set_ylabel("Noise (coefficient of variation squared)",fontsize=17)
	
	updateLogTicklabels(ax)	


def plotWriteLog(sheet, binned_x, binned_intrinsic, binned_extrinsic, binned_total, xlabel, directory): # plot and write mean expression vs. noise in logarithmic scale	
	fig = plt.figure(figsize = (6,4),dpi=300)
	ax = fig.add_subplot(111)

	xmin = float('inf')
	xmax = float('-inf')
	
	sheet.write(12,0,'Logarithmic scale')
	labels = ['Bin','# slices',' Normalized mean','2 std error -','2 std error +','Total noise','2 std error -','2 std error +',
		'Intrinsic noise','2 std error -','2 std error +','Extrinsic noise','2 std error -','2 std error +']
	for i in range(len(labels)):
		sheet.write(13,i,labels[i])
	
	mean_level = []
	total_noise = []
	intrinsic_noise = []
	extrinsic_noise = []
	background_log_tb = []
	background_log_im = []
	background_log_eg = []
	x_error_left = []
	x_error_right = []
	y_tot_error_up = []
	y_tot_error_down = []
	y_in_error_up = []
	y_in_error_down = []
	y_ex_error_up = []
	y_ex_error_down = []
	# Plot data
	for i in range(len(binned_x)):
		sheet.write(i+14,0,i+1)
		sheet.write(i+14,1,len(binned_x[i]))
		
		if len(binned_x[i])<=1:
			continue
		x = math.log10(numpy.mean(binned_x[i]))
		total = math.log10(numpy.mean(binned_total[i]))
		intrinsic = math.log10(numpy.mean(binned_intrinsic[i]))
		xtest=numpy.mean(binned_extrinsic[i]) # Some binned_extrinsic<0, we want to eliminate them
		if xtest<0:
			xtest=1	
		#extrinsic = math.log10(numpy.mean(binned_extrinsic[i]))
		extrinsic = math.log10(xtest)
		
		
		xe = numpy.std(binned_x[i])/math.sqrt(len(binned_x[i]))
		total_ye = numpy.std(binned_total[i])/math.sqrt(len(binned_total[i]))
		intrinsic_ye = numpy.std(binned_intrinsic[i])/math.sqrt(len(binned_intrinsic[i]))
		extrinsic_ye = numpy.std(binned_extrinsic[i])/math.sqrt(len(binned_extrinsic[i]))
		
		# Applying log to 2 standard errors
		# Transforming error bars in log scale: http://labs.physics.dur.ac.uk/skills/skills/logscales.php
		xe_left = abs(math.log10(abs(numpy.mean(binned_x[i])-2*xe))-x)
		xe_right = abs(math.log10(abs(numpy.mean(binned_x[i])+2*xe))-x)
		total_ye_lower = abs(math.log10(abs(numpy.mean(binned_total[i])-2*total_ye))-total)
		total_ye_upper = abs(math.log10(abs(numpy.mean(binned_total[i])+2*total_ye))-total)
		intrinsic_ye_lower = abs(math.log10(abs(numpy.mean(binned_intrinsic[i])-2*intrinsic_ye))-intrinsic)
		intrinsic_ye_upper = abs(math.log10(abs(numpy.mean(binned_intrinsic[i])+2*intrinsic_ye))-intrinsic)
		extrinsic_ye_lower = abs(math.log10(abs(numpy.mean(binned_extrinsic[i])-2*extrinsic_ye))-extrinsic)
		extrinsic_ye_upper = abs(math.log10(abs(numpy.mean(binned_extrinsic[i])+2*extrinsic_ye))-extrinsic)
			
		'''bg_log_tb = math.log10(background_bin_tb[i])
		bg_log_im = math.log10(background_bin_im[i])
		bg_log_eg = math.log10(background_bin_eg[i])'''
		
		xmin = min(xmin, x-xe_left)
		xmax = max(xmax, x+xe_right)

		mean_level.append(x)
		total_noise.append(total)
		intrinsic_noise.append(intrinsic)
		extrinsic_noise.append(extrinsic)
		'''background_log_tb.append(bg_log_tb)
		background_log_im.append(bg_log_im)
		background_log_eg.append(bg_log_eg)'''

		# Applying log to 2 standard errors
		# Transforming error bars in log scale: http://labs.physics.dur.ac.uk/skills/skills/logscales.php
		x_error_left.append(xe_left)
		x_error_right.append(xe_right)
		y_tot_error_down.append(total_ye_lower)
		y_tot_error_up.append(total_ye_upper)
		y_in_error_down.append(intrinsic_ye_lower)
		y_in_error_up.append(intrinsic_ye_upper)
		y_ex_error_down.append(extrinsic_ye_lower)
		y_ex_error_up.append(extrinsic_ye_upper)
		
		# Write current bin's data
		line = [x, xe_left, xe_right, total, total_ye_lower, total_ye_upper,
			intrinsic, intrinsic_ye_lower, intrinsic_ye_upper, extrinsic, extrinsic_ye_lower, extrinsic_ye_upper]
		for j in range(len(line)):
			sheet.write(i+14,j+2,line[j])
	
	# Plot				
	tot_pt = ax.scatter(mean_level, total_noise, s = 22, edgecolors='none', c=colors[0], label = "Total")
	in_pt = ax.scatter(mean_level, intrinsic_noise, s = 22, edgecolors='none', c=colors[1], label = "Intrinsic")
	ex_pt = ax.scatter(mean_level, extrinsic_noise, s = 22, edgecolors='none', c=colors[2], label = "Extrinsic")
	'''if (len(binned_x) == len(background_bin_tb)):
		bk_pt = ax.plot(mean_level, background_log_tb, ':', c = colors[0], label = "Baseline")
		ax.scatter(mean_level, background_log_tb, s = 22, edgecolors='none', c=colors[3])
		bk_pt = ax.plot(mean_level, background_log_im, ':', c = colors[1], label = "Baseline")
		ax.scatter(mean_level, background_log_im, s = 22, edgecolors='none', c=colors[3])
		bk_pt = ax.plot(mean_level, background_log_eg, ':', c = colors[2], label = "Baseline")
		ax.scatter(mean_level, background_log_eg, s = 22, edgecolors='none', c=colors[3])'''

	ax.errorbar(mean_level, total_noise, xerr = [x_error_left, x_error_right], yerr = [y_tot_error_down, y_tot_error_up], ls='none', c=colors[0],capsize=2, elinewidth=0.5)		
	ax.errorbar(mean_level, intrinsic_noise, xerr = [x_error_left,x_error_right], yerr = [y_in_error_down, y_in_error_up], ls='none', c=colors[1],capsize=2, elinewidth=0.5)
	ax.errorbar(mean_level, extrinsic_noise, xerr = [x_error_left,x_error_right], yerr = [y_ex_error_down, y_ex_error_up], ls='none', c=colors[2],capsize=2, elinewidth=0.5)
	
	# Legend
	bg_lg = mlines.Line2D([],[], color = colors[3], marker = 'o', markersize = 5, linestyle = '--')
	plt.legend([tot_pt,in_pt, ex_pt, bg_lg], ['Total', 'Intrinsic', 'Extrinsic', 'Baseline'], ncol = 4, \
	scatterpoints = 1, bbox_to_anchor = (0.5, 1.15), loc = 9, fontsize = 13)
	
	ax.tick_params(direction='in',labelsize = 14)
	xt = [1,2,4,6,8,10,20,40,60,80,100,200,400]
	for i in range(len(xt)):
		xt[i] = math.log10(xt[i])	
	ax.xaxis.set_ticks(xt)
	ax.set_xlim(findLogXlim(xmin,xmax)) # use minimum and maximum to find appropriate x-axis limits
	
	ylim = [-2, math.log10(0.80)] 
	yt = [0.01, 0.05, 0.10, 0.20, 0.40, 0.60]
	for i in range(len(yt)):
		yt[i] = math.log10(yt[i])
	ax.yaxis.set_ticks(yt)
	ax.set_ylim(ylim)	
	updateLogTicklabels(ax)	
	ax.set_ylabel("Noise",fontsize=22)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in grouped slices" % {'label':xlabel}, fontsize = 20)
	ax.tick_params(direction='in',labelsize = 14)
	# save figure
	fig.subplots_adjust(left=0.13, bottom=0.15, right=.95, top=.875, wspace=None, hspace=0.3)
	#fig.savefig(directory + "/logNoise_" + xlabel + ".tiff", format = "tiff", dpi=300)
	fig.savefig(directory + "/logNoise_" + xlabel + ".png", format = "png", dpi=300)

def plotWritesingle(sheet, binned_x, binned_total, xlabel, directory): # Qiyuan Add plot and write single gene mean expression vs. noise
	fig = plt.figure(figsize = (12,7),dpi=300)
	ax = fig.add_subplot(111)
	xmin = float('inf')
	xmax = float('-inf')
	ymax = float('-inf')
	
	mean_level = []
	total_noise = []
	x_error = []
	y_tot_error = []
	y_in_error = []
	y_ex_error = []
	
	# Calculate ys
	for i in range(len(binned_x)):
		xmin = min(xmin,numpy.mean(binned_x[i]))
		xmax = max(xmax,numpy.mean(binned_x[i])+2*numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])))
		ymax = max(ymax,numpy.mean(binned_total[i])+2*numpy.std(binned_total[i])/math.sqrt(len(binned_total[i])))
		mean_level.append(numpy.mean(binned_x[i]))
		total_noise.append(numpy.mean(binned_total[i]))
		x_error.append(2*numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])))
		y_tot_error.append(2*numpy.std(binned_total[i])/math.sqrt(len(binned_total[i])))
	# Plot
	'''if (len(binned_x) == len(background_bin_tb)):
		background_tb = ax.plot(mean_level, background_bin_tb, '--', c = colors[0], label="Baseline Noise")'''
	
	tot = ax.scatter(mean_level, total_noise, s = 22, edgecolors='none', c=colors[0], label = "Total")
	ax.legend(loc = 0, scatterpoints=1, fontsize = 14)
	'''if len(binned_x) == len(background_bin_tb):
		ax.scatter(mean_level, background_bin_tb, s = 22, edgecolors='none', c=colors[0])'''
	
	ax.errorbar(mean_level, total_noise, xerr = x_error, yerr = y_tot_error, ls='none', c=colors[0],capsize=2, elinewidth=0.5)
	
	ax.tick_params(direction='in',labelsize = 14)
	ax.xaxis.set_ticks(numpy.arange(0, xmax+10, determineTickInterval(xmax,10)))
	if xmax>=1: # need more space on the x-axis
		ax.set_xlim(0,xmax+xmin)
	else:
		ax.set_xlim(0,1)
	ax.yaxis.set_ticks(numpy.arange(0, ymax+.1, .2), determineTickInterval(ymax+ 0.1, 10))
	ax.set_ylim(0,ymax+0.05)
	ax.set_ylabel("Noise (coefficient of variation squared)", fontsize = 17)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in grouped slices" % {'label':xlabel}, fontsize = 10)

	# save figure
	fig.subplots_adjust(left=0.125, bottom=0.15, right=.95, top=.95, wspace=0.2, hspace=0.3)
	fig.savefig(directory + "/noise_" + xlabel + ".png", format = "png", dpi=300)
	# Write data
	labels = ['Bin','# slices','Normalized mean','Std error','Total noise','Std error']
	for i in range(len(labels)):
		sheet.write(0,i,labels[i])
	for i in range(len(binned_x)):
		line = [i+1, len(binned_x[i]), numpy.mean(binned_x[i]), numpy.std(binned_x[i])/math.sqrt(len(binned_x[i])),
			numpy.mean(binned_total[i]), numpy.std(binned_total[i])/math.sqrt(len(binned_total[i]))]
		for j in range(len(line)):
			sheet.write(i+1,j,line[j])

def plotRawHervsNoise(mean_level, total_noise, intrinsic_noise, extrinsic_noise, xlabel, directory):
	matplotlib.rcParams['xtick.minor.size'] = 0
	matplotlib.rcParams['xtick.minor.width'] = 0
	matplotlib.rcParams['ytick.minor.size'] = 0
	matplotlib.rcParams['ytick.minor.width'] = 0
	fig = plt.figure(figsize = (6,4),dpi=300)
	ax = fig.add_subplot(111)

	tot_pt = ax.scatter(mean_level, total_noise, s = 9, edgecolors = 'none', c=colors[0], label = "Total")
	#tot_patch = mpatches.Patch(color=colors[0], label="Total noise in single slices")
	ax.legend([tot_pt], ['Total noise in single slices'], ncol = 1, scatterpoints = 1, bbox_to_anchor = (0.5, 1.2), loc = 9, fontsize = 16)
	ax.set_xscale('log')
	ax.set_yscale('log')
    
	xmin, xmax = numpy.min(mean_level), numpy.max(mean_level)
	ax.set_xlim([0.9*xmin, 1.1*xmax])
	ymin, ymax = numpy.min(total_noise), numpy.max(total_noise)
	ax.set_ylim([0.9*ymin, 1.1*ymax])
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.set_yticks([0.01, 0.05, 0.2, 0.8, 2.0])
	ax.set_xticks([2,4,10,20,40,80,200])
	ax.tick_params(direction='in',labelsize = 14)
	ax.set_ylabel(r"Total Noise ($\textsf{CV}^{\textsf{\small{2}}}$)",fontsize=22)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in single slice" % {'label':xlabel}, fontsize = 20)
	fig.subplots_adjust(left=0.13, bottom=0.15, right=.95, top=.875, wspace=None, hspace=0.3)
	fig.savefig(directory + "/loglogHervsNoise_" + xlabel + ".png", format = "png", dpi=300)

def plotRawHer1vsNoise(mean_level, total_noise, intrinsic_noise, extrinsic_noise, xlabel, directory):
	matplotlib.rcParams['xtick.minor.size'] = 0
	matplotlib.rcParams['xtick.minor.width'] = 0
	matplotlib.rcParams['ytick.minor.size'] = 0
	matplotlib.rcParams['ytick.minor.width'] = 0
	fig = plt.figure(figsize = (6,4),dpi=300)
	ax = fig.add_subplot(111)

	tot_pt = ax.scatter(mean_level, total_noise, s = 9, edgecolors = 'none', c=colors[1], label = "Total")
	#tot_patch = mpatches.Patch(color=colors[0], label="Total noise in single slices")
	ax.legend([tot_pt], ['Total noise in single slices'], ncol = 1, scatterpoints = 1, bbox_to_anchor = (0.5, 1.2), loc = 9, fontsize = 16)
	ax.set_xscale('log')
	ax.set_yscale('log')
    
	xmin, xmax = numpy.min(mean_level), numpy.max(mean_level)
	ax.set_xlim([3*xmin, 1*xmax])
	ymin, ymax = numpy.min(mean_level), numpy.max(mean_level)
	ax.set_ylim([0.9*ymin, 0.04*ymax])
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.set_yticks([0.01, 0.05, 0.2, 0.8, 2.0])
	ax.set_xticks([2,4,10,20,40,80,200])
	ax.tick_params(direction='in',labelsize = 14)
	ax.set_ylabel(r"Total Noise ($\textsf{CV}^{\textsf{\small{2}}}$)",fontsize=22)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in single slice" % {'label':xlabel}, fontsize = 20)
	fig.subplots_adjust(left=0.13, bottom=0.15, right=.95, top=.875, wspace=None, hspace=0.3)
	fig.savefig(directory + "/loglogHer1vsNoise_" + xlabel + ".png", format = "png", dpi=300)

def plotRawHer7vsNoise(mean_level, total_noise, intrinsic_noise, extrinsic_noise, xlabel, directory):
	matplotlib.rcParams['xtick.minor.size'] = 0
	matplotlib.rcParams['xtick.minor.width'] = 0
	matplotlib.rcParams['ytick.minor.size'] = 0
	matplotlib.rcParams['ytick.minor.width'] = 0
	fig = plt.figure(figsize = (6,4),dpi=300)
	ax = fig.add_subplot(111)

	tot_pt = ax.scatter(mean_level, total_noise, s = 9, edgecolors = 'none', c=colors[2], label = "Total")
	#tot_patch = mpatches.Patch(color=colors[0], label="Total noise in single slices")
	ax.legend([tot_pt], ['Total noise in single slices'], ncol = 1, scatterpoints = 1, bbox_to_anchor = (0.5, 1.2), loc = 9, fontsize = 16)
	ax.set_xscale('log')
	ax.set_yscale('log')
    
	xmin, xmax = numpy.min(mean_level), numpy.max(mean_level)
	ax.set_xlim([0.9*xmin, 1.1*xmax])
	ymin, ymax = numpy.min(mean_level), numpy.max(mean_level)
	ax.set_ylim([0.9*ymin, 0.04*ymax])
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.set_yticks([0.01, 0.05, 0.2, 0.8, 2.0])
	ax.set_xticks([2,4,10,20,40,80,200])
	ax.tick_params(direction='in',labelsize = 14)
	ax.set_ylabel(r"Total Noise ($\textsf{CV}^{\textsf{\small{2}}}$)",fontsize=22)
	ax.set_xlabel(r"Mean $\textit{%(label)s}$ mRNA in single slice" % {'label':xlabel}, fontsize = 20)
	fig.subplots_adjust(left=0.13, bottom=0.15, right=.95, top=.875, wspace=None, hspace=0.3)
	fig.savefig(directory + "/loglogHer1vsNoise_" + xlabel + ".png", format = "png", dpi=300)

def writeSPSS(ws, slice_mean_her, slice_intrinsic_noise, slice_extrinsic_noise): # write data for SPSS statistical analysis
	labels = ['Low expression noise type','Noise level','','High expression noise type','Noise level']	
	for i in range(len(labels)):
		ws.write(0,i,labels[i])	
	threshold = (max(slice_mean_her)-min(slice_mean_her))/2
	low_index = 1
	high_index = 1	
	for i in range(len(slice_mean_her)):	
		if slice_mean_her[i]<threshold:
			ws.write(low_index,0,1)
			ws.write(low_index,1,slice_intrinsic_noise[i])
			ws.write(low_index+1,0,2)
			ws.write(low_index+1,1,slice_extrinsic_noise[i])
			low_index+=2
		else:
			ws.write(high_index,3,1)
			ws.write(high_index,4,slice_intrinsic_noise[i])
			ws.write(high_index+1,3,2)
			ws.write(high_index+1,4,slice_extrinsic_noise[i])
			high_index+=2				

	
def write_aggregate_data(directory, slice_mean_her, slice_mean_her1, slice_mean_her7, slice_mean_hm, slice_intrinsic_noise, \
slice_extrinsic_noise, slice_total_noise, all_her1, all_her7, num_bins):
	workbook = xlwt.Workbook(encoding="ascii")
	# 0. Create workshee "aggregate_data"
	ws = workbook.add_sheet("aggregate_data")
	# 1. Write headers
	headers = ["avg_her", "avg_her1", "avg_her7", "avg_harmonic_mean", "intrinsic", "extrinsic", "total_noise",\
	"num_cells", "all her1...", "all her7..."]
	for i in range(len(headers)):
		ws.write(0, i, headers[i])
	# 2. Write all rows of data, each row is a slice
	for i in range(len(slice_mean_her)):
		ws.write(i + 1, 0, slice_mean_her[i]) # avg_her
		ws.write(i + 1, 1, slice_mean_her1[i]) # avg_her1
		ws.write(i + 1, 2, slice_mean_her7[i]) # avg_her7
		ws.write(i + 1, 3, slice_mean_hm[i]) # avg_harmonic_mean
		ws.write(i + 1, 4, slice_intrinsic_noise[i]) # intrinsic noise levels of slices
		ws.write(i + 1, 5, slice_extrinsic_noise[i]) # extrinsic noise levels of slices
		ws.write(i + 1, 6, slice_total_noise[i]) # total noise levels of slices
		ws.write(i + 1, 7, len(all_her1[i])) # number of cells in this slice
		for j in range(len(all_her1[i])): # for each cell in this slice, wirte all her1 first, then all her7
			ws.write(i + 1, 8 + j, all_her1[i][j])
			ws.write(i + 1, 8 + j + len(all_her1[i]), all_her7[i][j])
	# 3. Create a new worksheet called "bin_bounds".
	# Write into that sheet the upper bound of total her levels in each bin
	bound_ws = workbook.add_sheet("bin_bounds")
	min_her = min(slice_mean_her)
	max_her = max(slice_mean_her)
	her_interval = (max_her - min_her) / num_bins
	for i in range(num_bins):
		bound_ws.write(0, i, "Bin " + str(i) + "her_upper bound")
		bound_ws.write(1, i , min_her + (i + 1) * her_interval)
	workbook.save(directory + "/combined_slices.xls")

def write_histogram_input(directory, slice_mean_her, all_her1, all_her7):
	workbook = xlwt.Workbook(encoding="ascii")
	
	# 0. Create worksheet "histogram"
	ws = workbook.add_sheet("histogram")
	
	# 1. Write headers
	headers = ["Her 1","Her 7","Total Her"]
	for i in range(len(headers)):
		ws.write(0, i, headers[i])
		
	# 2. Write all rows of data, each row is a slice
	curR = 1
	for i in range(len(slice_mean_her)):
		for j in range(len(all_her1[i])): # for each cell in this slice, wirte all her1 first, then all her7
			ws.write(curR, 0, all_her1[i][j])
			ws.write(curR, 1, all_her7[i][j])
			ws.write(curR, 2, all_her1[i][j] + all_her7[i][j])
			curR = curR + 1
	workbook.save(directory + "/histogram_her.xls")
	
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print('plot_noise.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print('plot_noise.py: Number of embryos must be larger than zero.')
		exit(1)
	if not shared.isInt(sys.argv[2]):
		print('plot_noise.py: Number of bins must be an integer.')
		exit(1)
	elif int(sys.argv[2])<=0:
		print('plot_noise.py: Number of bins must be larger than zero.')
		exit(1)
	
	num_embryos = int(sys.argv[1])
	num_bins = int(sys.argv[2])	# fix bin size of 5 is defined at the beginning in this script
        
	if len(sys.argv)==num_embryos+4:
		inputs = sys.argv[3:num_embryos+3]
		directory = sys.argv[num_embryos+3]		
	else:
		usage()
	
	shared.ensureDir(directory)
	
	# mean levels of all slices of all regions of all embryos [slice]
	slice_mean_her1 = []
	slice_mean_her7 = []	
	slice_mean_her = []
	slice_mean_hm = [] # harmonic mean
	slice_intrinsic_noise = []
	slice_extrinsic_noise = []
	slice_total_noise = []
	slice_her1_total_noise = []
	slice_her7_total_noise = []
	# 2D: [slice][cell]: matrix to store the levels of her1 and her7 of all cells of all slices in all regions in all embryos. 
	# The slice index matches slice_mean_... things above
	all_her1 = []
	all_her7 = []
		
	for i in range(len(inputs)): # embryo
		# Open embryo data
		if not os.path.isfile(inputs[i]):
			print('plot_noise.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print('plot_noise.py: Cannot open file "'+inputs[i]+'".')
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
					for l in range(num_cells): # background subtraction
						# Take the cell's data only if its expression levels are positive after background subtraction	
						if row[8+2*l].value>0 and row[8+2*l+1].value>0:
							her1.append(row[8+2*l].value)
							her7.append(row[8+2*l+1].value)
					
					if len(her1)>=3: # take data only if 3 or more cells
						her1_mean = numpy.mean(her1)
						her7_mean = numpy.mean(her7)
						# Put data into structures for all slices of all regions of all embryos												
						slice_mean_her1.append(her1_mean) # her1 mean
						slice_mean_her7.append(her7_mean) # her7 mean
						slice_mean_her.append(her1_mean+her7_mean) # her mean (her1 + her7)
						slice_mean_hm.append(2/(1/her1_mean+1/her7_mean)) # harmonic mean of her1 and her7
						all_her1.append(her1)
						all_her7.append(her7)
						
						# Intrinsic noise
						intrinsic_noise = 0						
						for m in range(len(her1)):
							intrinsic_noise += (her1[m]/her1_mean - her7[m]/her7_mean)**2
						intrinsic_noise = intrinsic_noise / len(her1) / 2
						slice_intrinsic_noise.append(intrinsic_noise)
						
						# Extrinsic noise
						product = []
						for l in range(len(her1)):
							product.append(her1[l]*her7[l])
						extrinsic_noise = (numpy.mean(product) - her1_mean*her7_mean) / (her1_mean*her7_mean)
						slice_extrinsic_noise.append(extrinsic_noise)
						
						# Total noise (sum of intrinsic and extrinsic noises)
						slice_total_noise.append(intrinsic_noise+extrinsic_noise)

						# Qiyuan Add: Total noise for each gene
						her1_total_noise = 0 # her1 single gene total noise CV2
						for m in range(len(her1)):
							her1_total_noise = (numpy.std(her1)/her1_mean)**2
						slice_her1_total_noise.append(her1_total_noise)
						
						her7_total_noise = 0 # her7 single gene total noise CV2
						for m in range(len(her7)):
							her7_total_noise = (numpy.std(her7)/her7_mean)**2
						slice_her7_total_noise.append(her7_total_noise)

	# Plot raw mean_her vs noise before binned
	plotRawHervsNoise(slice_mean_her, slice_total_noise, slice_intrinsic_noise, slice_extrinsic_noise, "her", directory)

	# Qiyuan Add Plot raw mean_her1 vs her1 single gene total noise before binned
	plotRawHer1vsNoise(slice_mean_her1, slice_her1_total_noise, slice_intrinsic_noise, slice_extrinsic_noise, "her1", directory)
	
	# Qiyuan Add Plot raw mean_her7 vs her7 single gene total noise before binned
	plotRawHer7vsNoise(slice_mean_her7, slice_her7_total_noise, slice_intrinsic_noise, slice_extrinsic_noise, "her7", directory)

	# Bin data				
	binned_data_her1 = binData_fix(slice_mean_her1, slice_intrinsic_noise, slice_extrinsic_noise, slice_total_noise)			
	binned_data_her7 = binData_fix(slice_mean_her7, slice_intrinsic_noise, slice_extrinsic_noise, slice_total_noise)	
	binned_data_her = binData_fix(slice_mean_her, slice_intrinsic_noise, slice_extrinsic_noise, slice_total_noise)	
	binned_data_hm = binData_fix(slice_mean_hm, slice_intrinsic_noise, slice_extrinsic_noise, slice_total_noise)
	binned_single_her1_data = binData_single_gene_fix(slice_mean_her1,slice_her1_total_noise) #Qiyuan Add
	binned_single_her7_data = binData_single_gene_fix(slice_mean_her7,slice_her7_total_noise) #Qiyuan Add

	# Write data of all slices and cells' her1 and her7 levels
	write_aggregate_data(directory, slice_mean_her, slice_mean_her1, slice_mean_her7, slice_mean_hm, slice_intrinsic_noise, \
	slice_extrinsic_noise, slice_total_noise, all_her1, all_her7, num_bins)
	#write_histogram_input(directory, slice_mean_her, all_her1, all_her7)
	
	# Plot and write data of binned noise and concentration
	workbook = xlwt.Workbook(encoding="ascii")
	her1_worksheet = workbook.add_sheet("Her1")
	plotWrite(her1_worksheet, binned_data_her1[0], binned_data_her1[1], binned_data_her1[2], binned_data_her1[3], "her1", directory)
	her7_worksheet = workbook.add_sheet("Her7")
	plotWrite(her7_worksheet, binned_data_her7[0], binned_data_her7[1], binned_data_her7[2], binned_data_her7[3], "her7", directory)
	her_worksheet = workbook.add_sheet("Her")
	plotWrite(her_worksheet, binned_data_her[0], binned_data_her[1], binned_data_her[2], binned_data_her[3], "her", directory)
	hm_worksheet = workbook.add_sheet("Harmonic mean")
	plotWrite(hm_worksheet, binned_data_hm[0], binned_data_hm[1], binned_data_hm[2], binned_data_hm[3], "Harmonic mean", directory)
	her1singlenoise_worksheet = workbook.add_sheet("Her1 singl noise")
	plotWritesingle(her1singlenoise_worksheet, binned_single_her1_data[0],binned_single_her1_data[1],"Her1 singl noise", directory)
	her7singlenoise_worksheet = workbook.add_sheet("Her7 singl noise")
	plotWritesingle(her7singlenoise_worksheet, binned_single_her7_data[0],binned_single_her7_data[1],"Her7 singl noise", directory)

	plotWriteLog(her1_worksheet, binned_data_her1[0], binned_data_her1[1], binned_data_her1[2], binned_data_her1[3], "her1", directory)
	plotWriteLog(her7_worksheet, binned_data_her7[0], binned_data_her7[1], binned_data_her7[2], binned_data_her7[3], "her7", directory)
	plotWriteLog(her_worksheet, binned_data_her[0], binned_data_her[1], binned_data_her[2], binned_data_her[3], "her", directory)
	plotWriteLog(hm_worksheet, binned_data_hm[0], binned_data_hm[1], binned_data_hm[2], binned_data_hm[3], "Harmonic Mean", directory)
	
	# Write data for spss	
	spss_worksheet = workbook.add_sheet("spss_onewayANOVA")	
	writeSPSS(spss_worksheet, slice_mean_her, slice_intrinsic_noise, slice_extrinsic_noise)
	
	binspss_worksheet = workbook.add_sheet("spss_bin")	
	writeSPSS5_bin(binspss_worksheet, binned_data_her[0], binned_data_her[1], binned_data_her[2], binned_data_her[3],len(bin_range))
	
	binher1spss_worksheet = workbook.add_sheet("spss_bin_her1_single")	
	writeSPSS5_single_bin(binher1spss_worksheet, binned_single_her1_data[0],binned_single_her1_data[1],len(bin_range)) #Oriana add

	binher7spss_worksheet = workbook.add_sheet("spss_bin_her7_single")	
	writeSPSS5_single_bin(binher7spss_worksheet, binned_single_her7_data[0],binned_single_her7_data[1],len(bin_range)) # Oriana add



	workbook.save(directory + "/noise.xls")
		
def usage():
	print("plot_noise.py: Invalid command-line arguments")
	print("Format: python plot_noise.py <number of embryos> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print("Example: python plot_noise.py 20 ../wildtypefulldataset/output/embryo1/slices.xls \
	../wildtypefulldataset/output/embryo2/slices.xls ... ../wildtypefulldataset/output/embryo20/slices.xls \
	../wildtypefulldataset/output")
	exit(1)

main()
