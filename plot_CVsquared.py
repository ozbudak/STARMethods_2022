"""
Plot spatial noise: cell position vs. CV^2 
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

gene_colors = ['b','r','g'] # her1, her7, her
expression_colors = ['#2A6EFF','m','#008d74'] # low (sky blue), medium (pink/purple), high (green)
region_colors = ['#2A6EFF','m','#008d74'] # posterior (sky blue), middle (pink/purple), anterior (green)

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2,5,10,20,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def updateXTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()] # int
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=5)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()] # float
	ax.set_yticklabels(ylabels)

def separateThreeGroups(all_mean, slice_mean, slice_cv_squared, num_slices): # separate data into three groups based on expression level
	thresholds = []	
	groups_cv_squared = [] # three-dimensional array 
	groups_mean = [] # three-dimensional array
	
	# Get the threshold for each group = average of 3 consecutive increasing mean levels of mRNAs
	num_slices_in_group = int(math.floor(len(all_mean)/3))		
	all_mean = sorted(all_mean)
	for i in range(2):
		thresholds.append(numpy.mean(all_mean[num_slices_in_group*(i+1):num_slices_in_group*(i+1)+2]))
	thresholds.append(max(all_mean))
		
	# Set up three-dimensional arrays [slice][group index: low, medium, high][embryo in that slice having mRNA levels in that group]
	for i in range(num_slices):	
		groups_cv_squared.append([])	
		groups_mean.append([])	
		for j in range(3):
			groups_cv_squared[i].append([])			
			groups_mean[i].append([])
	
	# Separate slices into three groups
	for i in range(num_slices):
		for j in range(len(slice_mean[i])):
			for k in range(len(thresholds)):		
				if slice_mean[i][j] <= thresholds[k]:
					groups_cv_squared[i][k].append(slice_cv_squared[i][j])
					groups_mean[i][k].append(slice_mean[i][j])
					break
					
	return groups_cv_squared, groups_mean			

def plot_grouped(ax, slice_cv_squared, groups_cv_squared, num_slices): # plot cell position vs. CV^2 with three separate expression groups
	ymax = 0
	for j in range(num_slices): # cell position		
		for k in range(3): # low, medium, high expression
			if len(groups_cv_squared[j][k])>0: # enough data to plot					
				ax.scatter(j, numpy.mean(groups_cv_squared[j][k]), c=expression_colors[k], s = 22, edgecolors='none')
				ax.errorbar(j, numpy.mean(groups_cv_squared[j][k]),
					yerr=2*numpy.std(groups_cv_squared[j][k])/math.sqrt(len(groups_cv_squared[j][k])), ls='none',c=expression_colors[k], capsize=2)
				ymax = max(ymax, numpy.mean(groups_cv_squared[j][k])+2*numpy.std(groups_cv_squared[j][k])/math.sqrt(len(groups_cv_squared[j][k])))
			
	ax.xaxis.set_ticks(numpy.arange(0,num_slices+1,determineTickInterval(num_slices,6)))
	ax.set_xlim(-1,num_slices)
	ax.set_xlabel(r"Cell position (posterior - anterior)")
	ax.yaxis.set_ticks(numpy.arange(0,ymax+1,determineTickInterval(ymax+1,6)))
	ax.set_ylim(0,ymax*1.05)
	ax.tick_params(direction='in')
	updateXTicklabels(ax)			
		
	ax.set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her)}")
	
	# Figure legend
	low = plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=expression_colors[0], \
	markeredgecolor=expression_colors[0], markersize=4)
	medium = plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=expression_colors[1], \
	markeredgecolor=expression_colors[1], markersize=4)
	high = plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=expression_colors[2], \
	markeredgecolor=expression_colors[2], markersize=4)	
	ax.legend([low,medium,high],[r"Low expression",r"Medium expression",r"High expression"], \
	numpoints=1, loc=9, bbox_to_anchor=(0.45,1.15), ncol=3, fontsize = 10, \
	handletextpad = 0.05, labelspacing = 0.05, columnspacing = 0.08)

def plot(ax, gene, slice_cv_squared,num_slices): # plot cell position vs. CV^2 
	for j in range(num_slices):
		ax.scatter(j, numpy.mean(slice_cv_squared[gene][j]), c='#722AFF', s = 22, edgecolors='none')
		ax.errorbar(j, numpy.mean(slice_cv_squared[gene][j]),
			yerr=2*numpy.std(slice_cv_squared[gene][j])/math.sqrt(len(slice_cv_squared[gene][j])), ls='none', c='#722AFF',capsize=2)
	ax.xaxis.set_ticks(numpy.arange(0,num_slices+1,determineTickInterval(num_slices,6)))
	ax.set_xlim(-1,num_slices)
	ax.set_xlabel(r"Cell position (posterior - anterior)")
	start,end = ax.get_ylim()
	ax.yaxis.set_ticks(numpy.arange(0,end+1,determineTickInterval(end,5)))
	ax.set_ylim(0,end)
	ax.tick_params(direction='in')
	updateXTicklabels(ax)

def write(sheet, slice_mean, slice_cv_squared, groups_mean, groups_cv_squared, num_slices, num_groups): # write data to Excel
	labels = ['Cell position','Mean','CV^2','Std error','# slices','', 
		'Low mean','Std error','CV^2 mean','Std error','# slices','',
		'Medium mean','Std error','CV^2 mean','Std error','# slices','',
		'High mean','Std error','CV^2 mean','Std error','# slices','',
		'Raw data (CV^2)']
	for i in range(len(labels)):
		sheet.write(0,i,labels[i])
	
	for i in range(num_slices): # cell position
		line = []
		if len(slice_cv_squared[i])>0:
			line = [i, numpy.mean(slice_mean[i]), numpy.mean(slice_cv_squared[i]), numpy.std(slice_cv_squared[i])/math.sqrt(len(slice_cv_squared[i])),len(slice_mean[i]),""]
			for j in range(num_groups): # add data based on three expression groups to the right
				if len(groups_mean[i][j])>0: # enough data to write
					line.append(numpy.mean(groups_mean[i][j]))
					line.append(numpy.std(groups_mean[i][j])/math.sqrt(len(groups_mean[i][j])))
					line.append(numpy.mean(groups_cv_squared[i][j]))
					line.append(numpy.std(groups_cv_squared[i][j])/math.sqrt(len(groups_cv_squared[i][j])))
					line.append(len(groups_mean[i][j]))
					line.append("")
				else:
					line += [""]*4 + [0] + [""]

		for j in range(len(line)):
			sheet.write(i+1,j,line[j])
		# we comment out the following piece of code because this raw data is rather unnecessary. 
		# we will write raw data for her mRNA levels into separate file
		'''
		for j in range(len(slice_cv_squared[i])):
			sheet.write(i+1,j+6+num_groups*6,slice_cv_squared[i][j]) # write all CV^2 values within the same cell position
		'''
		
def writeSPSS(spss_worksheets,slice_cv_squared): # write data to Excel for SPSS statistical analysis
	labels = ['Region(L/R)','CV^2']
	for m in range(3): # her1/her7/her
		ws = spss_worksheets[m]
		for (i, label) in enumerate(labels):
			ws.write(0, i, label)
		row_index = 1
		for i in range(len(slice_cv_squared[m])): # cell position
			for j in range(len(slice_cv_squared[m][i])):			
				# Posterior vs. anterior
				ws.write(row_index,0,1) if i<len(slice_cv_squared[m])/2 else ws.write(row_index,0,2)
				ws.write(row_index,1,slice_cv_squared[m][i][j]) 
				row_index += 1

def plot_heatmap_cv_her_pos(slice_cv2, slice_mean, state_name, save_file_name, num_slices):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	heat_map_colors = ['#FF7A33', '#FFF333', '#C1FF33']
	current_pos_bin = 0 
	for i in range(num_slices):
		if (i >= int(float(num_slices) * float((current_pos_bin + 1)) / float(3))):
			current_pos_bin += 1
		if len(slice_cv2[i]) > 0:
			ax.scatter(slice_mean[i], slice_cv2[i], color = heat_map_colors[current_pos_bin], alpha = 0.2, marker = '.')
	ax.set_xlim(0, 200)
	ax.set_xlabel('$\mathit{' + state_name + '}$ mRNAs')
	ax.set_ylabel('$CV^2$') 
	fig.savefig(save_file_name, format = "png", dpi = 300)

def plot_heatmap_binned_cvHerPos(slice_cv2, slice_mean, state_name, save_file_name, \
num_slices, interval, max_her, write_spss_ws, write_sum_ws):
	# create new excel file
	labels = ["bin_index", "Pos=0/Ant=2", "cv2", "mean"]
	for i, label in enumerate(labels):
		write_spss_ws.write(0,i,label)
	labels = ["bin_index", "Pos=0/Ant=2", "mean_her", "mean_cv2", "her_ste", "cv2_ste"]
	for i, label in enumerate(labels):
		write_sum_ws.write(0,i,label)
		
	# create new plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	heatmap_regions = 3
	current_pos_bin = 0 
	num_bins = int(max_her / interval)
	binned_mean = [] #[bin][region][data point]
	binned_cv2 = [] 
	for i in range(num_bins):
		binned_mean.append([[], [], []]) # posterior, middle, anterior
		binned_cv2.append([[], [], []])
	
	### BIN DATA 
	for i in range (num_slices):# slice
		# if one third of the embryo is over, we move onto a new region
		if (i >= int(float(num_slices) * float((current_pos_bin + 1)) / float(3))):
			current_pos_bin += 1
		for j in range(len(slice_cv2[i])):# embryo having such a slice
			bin_index = 0
			add = True
			while bin_index < num_bins:
				if (slice_mean[i][j] <= ((bin_index + 1) * interval)):
					break
				else: 
					bin_index += 1
				if bin_index ==  num_bins:
					add = False
			if add: 
				(binned_mean[bin_index][current_pos_bin]).append(slice_mean[i][j])
				(binned_cv2[bin_index][current_pos_bin]).append(slice_cv2[i][j])
	
	### WRITE spss DATA
	row_index = 1
	for i in range(num_bins):
		for j in range(heatmap_regions):
			if (len(binned_mean[i][j]) > 0) and j != 1: # right now we only write down values of posterior and anterior, we don't need the middle part
				for k in range(len(binned_cv2[i][j])):
					line = [i, j, binned_cv2[i][j][k], binned_mean[i][j][k]]
					for l, data_item in enumerate(line):
						write_spss_ws.write(row_index, l, data_item)
					row_index += 1
						
	### PLOT and WRITE plot DATA
	xmax = -float('inf')
	ymax = -float('inf')
	row_index = 1
	for i in range(num_bins):
		for j in range(heatmap_regions):
			if (len(binned_mean[i][j]) > 0) and j != 1: # right now we skip plotting the middle part of the embryo, we only plot the posterior and the anterior
				# calculate number
				mean = numpy.mean(binned_mean[i][j])
				cv2 = numpy.mean(binned_cv2[i][j])
				xerr = numpy.std(binned_mean[i][j]) / math.sqrt(len(binned_mean[i][j])) # one standard error
				yerr = numpy.std(binned_cv2[i][j]) / math.sqrt(len(binned_cv2[i][j]))
				# plot
				ax.scatter(mean, cv2, color = region_colors[j])
				ax.errorbar(mean, cv2, xerr = 2 * xerr, \
				yerr = 2 * yerr, color = region_colors[j], capsize=2)
				xmax = max(xmax, mean + xerr)
				ymax = max(ymax, cv2 + yerr)
				# write summary data
				line = [i, j, mean, cv2, xerr, yerr]
				for k, line_data in enumerate(line):
					write_sum_ws.write(row_index, k, line_data)
				row_index += 1
					
	### Fix the figures' ticks
	ax.set_xlim (0, xmax + 5)
	ax.set_ylim (0, ymax + 0.1)
	ax.set_xticks([(i * interval) for i in range(num_bins + 1)]) 
	ax.set_yticks(numpy.arange(0, ymax + 0.1, determineTickInterval(ymax, 5)))
	updateXTicklabels(ax)
	ax.set_xlabel('Total $\mathit{' + state_name + '}$ mRNA')
	ax.set_ylabel('$CV^2$') 
	# Figure legend
	posterior = plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=region_colors[0], markeredgecolor=region_colors[0], markersize=8)
	anterior = plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=region_colors[2], markeredgecolor=region_colors[2], markersize=8)	
	ax.legend([posterior, anterior], [r"Posterior",r"Anterior"], loc=9, bbox_to_anchor=(0.5,1.15), ncol=2, numpoints=1, fontsize=12)
	# adjust the plot
	fig.subplots_adjust(left=0.1, bottom=0.1, right=.975, top=.85,  wspace=None, hspace=None)
	fig.savefig(save_file_name, format = "png", dpi = 300)

def write_raw_her_cv2(wb, groups_mean_her, group_cv_squared_her, num_slices, num_groups, group_names):
	assert num_groups == len(group_names), "plot_CVsquared.py: Number of groups is not the same as the number of sheets names provided"
	labels = ["slice_index", "cv2"]

	for i in range(num_groups):
		row_index = 1
		# create worksheet and write labels
		ws = wb.add_sheet(group_names[i])
		for j, label in enumerate(labels):
			ws.write(0, j, label)
		# write data
		for k in range(num_slices): # for each slide
			cv2_in_group = group_cv_squared_her[k][i] # get the data of this group in this slide
			for j, cv2 in enumerate(cv2_in_group):
				ws.write(row_index, 0, k) # write slide index
				ws.write(row_index, 1, cv2) # write cv2 in this group
				row_index += 1
		
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print('plot_CVsquared.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print('plot_CVsquared.py: Number of embryos must be larger than zero.')
		exit(1)		
	
	num_embryos = int(sys.argv[1])
	    	        
	if len(sys.argv)==num_embryos+3:
		inputs = sys.argv[2:num_embryos+2]
		directory = sys.argv[num_embryos+2]		
	else:
		usage()
			
	shared.ensureDir(directory)
		       
	slice_mean_her1 = [] # two-dimensional array [slice][embryo having that slice]
	slice_mean_her7 = []	
	slice_mean_her = []	
	slice_cv_squared_her1 = [] # two-dimensional array [slice][embryo having that slice]
	slice_cv_squared_her7 = []
	slice_cv_squared_her = []
	
	num_slices = []	# 1D: Store all mean mRNA levels of all slices, regradless of the position
	all_mean_her1 = []
	all_mean_her7 = []
	all_mean_her = []
		
	for i in range(len(inputs)): # embryo
		# Open embryo data
		if not os.path.isfile(inputs[i]):
			print('plot_CVsquared.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print('plot_CVsquared.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		worksheets = workbook.sheets()	
		for j in range(len(worksheets)): # region 
			worksheet = worksheets[j]
			file_len =  worksheet.nrows		
			num_slices.append(file_len-1)	
						
			for k in range(1,file_len): # slice (cell position)
				row = list(worksheet.row(k))				
				
				# Set up multidimensional array
				if len(slice_mean_her1)<k:
					slice_mean_her1.append([])
					slice_cv_squared_her1.append([])
				if len(slice_mean_her7)<k:
					slice_mean_her7.append([])
					slice_cv_squared_her7.append([])
				if len(slice_mean_her)<k:
					slice_mean_her.append([])
					slice_cv_squared_her.append([])	
					
				if isinstance(row[1].value,float): # valid slice
					num_cells = int(row[1].value) # number of cells within this slice						
					her1 = []
					her7 = []
					her = []
					
					for l in range(num_cells): # background subtraction					
						# Take the cell's data only if its expression levels are positive after background subtraction	
						if row[8+2*l].value>0 and row[8+2*l+1].value>0:
							her1.append(row[8+2*l].value)
							her7.append(row[8+2*l+1].value)
							her.append(row[8+2*l].value + row[8+2*l+1].value)
					
					if len(her1)>=3: # valid only if there are more than 2 cells						
						# Mean
						her1_mean = numpy.mean(her1)
						her7_mean = numpy.mean(her7)
						her_mean = numpy.mean(her)						
						slice_mean_her1[k-1].append(her1_mean) # store mean for this cell position
						slice_mean_her7[k-1].append(her7_mean)
						slice_mean_her[k-1].append(her_mean)						
						all_mean_her1.append(her1_mean) # store mean regardless of cell position
						all_mean_her7.append(her7_mean)
						all_mean_her.append(her_mean)
						
						# CV^2 
						slice_cv_squared_her1[k-1].append((numpy.std(her1)/her1_mean)**2)
						slice_cv_squared_her7[k-1].append((numpy.std(her7)/her7_mean)**2)
						slice_cv_squared_her[k-1].append((numpy.std(her)/her_mean)**2)

	# Determine number of slices with at least 80% of embryos for analysis
	nSlices = sorted(num_slices)[int(num_embryos*0.2)]
	
	# Divide data into three groups based on average RNA level ---> 3D arrays [slice][group][embryo]
	groups_cv_squared_her1,groups_mean_her1 = separateThreeGroups(all_mean_her1,slice_mean_her1,slice_cv_squared_her1,nSlices)
	groups_cv_squared_her7,groups_mean_her7 = separateThreeGroups(all_mean_her7,slice_mean_her7,slice_cv_squared_her7,nSlices)
	groups_cv_squared_her,groups_mean_her = separateThreeGroups(all_mean_her,slice_mean_her,slice_cv_squared_her,nSlices)
	
	slice_cv_squared = [slice_cv_squared_her1,slice_cv_squared_her7,slice_cv_squared_her]
	groups_cv_squared = [groups_cv_squared_her1,groups_cv_squared_her7,groups_cv_squared_her]	

	# Plot cell position vs. CV^2 with three separate expression groups (low, medium, high)
	fig = plt.figure(figsize=(4.5,4),dpi=300)
	plot_grouped(fig.add_subplot(111),slice_cv_squared_her, groups_cv_squared_her, nSlices)
	fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.9,  wspace=None, hspace=0.25)	
	fig.savefig(directory + "/" + "CVsquared_grouped_her.png", format = "png", dpi=300)
	
	# Plot cell position vs. CV^2 with all expression groups combined in her1
	fig = plt.figure(figsize=(4.5,4),dpi=300)
	ax_her1 = fig.add_subplot(111)
	plot(ax_her1, 0, slice_cv_squared, nSlices)	
	ax_her1.set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her1)}")	
	fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.9,  wspace=None, hspace=.3)	
	fig.savefig(directory + "/" + "CVsquared_her1.png", format = "png", dpi=300)
	
	# Plot cell position vs. CV^2 with all expression groups combined in her7	
	fig = plt.figure(figsize=(4.5,4),dpi=300)
	ax_her7 = fig.add_subplot(111)
	plot(ax_her7, 1, slice_cv_squared, nSlices)	
	ax_her7.set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ \textit{(her7)}")	
	fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.9,  wspace=None, hspace=.3)	
	fig.savefig(directory + "/" + "CVsquared_her7.png", format = "png", dpi=300)

	# Write data to Excel
	workbook = xlwt.Workbook(encoding="ascii")
	write(workbook.add_sheet("Her1"), slice_mean_her1, slice_cv_squared_her1, groups_mean_her1, groups_cv_squared_her1, nSlices, 3)
	write(workbook.add_sheet("Her7"), slice_mean_her7, slice_cv_squared_her7, groups_mean_her7, groups_cv_squared_her7, nSlices, 3)
	write(workbook.add_sheet("Her"), slice_mean_her, slice_cv_squared_her, groups_mean_her, groups_cv_squared_her, nSlices, 3)
	write_raw_her_cv2(workbook, groups_mean_her, groups_cv_squared_her, nSlices, 3, ["her_low", "her_medium", "her_high"])

	# Now, slice_cv_squared_her and slice_mean_her and their lot will be 2D arrays, in the following format: 
	# [slice][embryo has this slice]
	# create binned_cv_her_heatmap.png
	interval = 15
	max_her = 120
	plot_heatmap_binned_cvHerPos(slice_cv_squared_her, slice_mean_her, 'her', \
	directory + '/binned_cv_her_heatmap.png', nSlices, interval, \
	max_her, workbook.add_sheet("binned_cv2_her_pos_spss"), workbook.add_sheet("binned_cv2_her_pos_summary"))
	
	# Write data for SPSS statistical analysis
	spss_worksheets = [workbook.add_sheet("spss_ANOVA_her1"), workbook.add_sheet("spss_ANOVA_her7"), workbook.add_sheet("spss_ANOVA_her")]
	writeSPSS(spss_worksheets,slice_cv_squared) 	
	workbook.save(directory + "/CVsquared.xls")
		
def usage():
	print("plot_CVsquared.py: Invalid command-line arguments")
	print("Format: python plot_CVsquared.py <number of embryos> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print("Example: python plot_CVsquared.py 20 wildtypefulldataset/embryo1/slices.xls wildtypefulldataset/embryo2/slices.xls ... wildtypefulldataset/embryo20/slices.xls wildtypefulldataset")
	exit(1)

main()
