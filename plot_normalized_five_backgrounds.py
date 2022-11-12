
import sys, shared, os
import numpy, math
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as mpatches
import xlrd, xlwt
from xlrd import XLRDError
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering
width = 0.5 # bar width for plotting

## index of the columns that contain information we need from input file (raw_noise.xls)
tot_input_index = 7
in_input_index = 5
ex_input_index = 6

# indices of all kinds of noises, used for common data structure
intrinsic = 0 
extrinsic = 1
total = 2

# Name of noises
noise_name = ["Intrinsic", "Extrinsic", "Total"]
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
	
def calculate_noise_data(num_backgrounds, inputs, ws_sum, ws_spss):
	#write labels into the excel file
	labels = ["Genetic background", "raw_in_noise", "raw_ex_noise", "raw_tot_noise", "normalized_in_noise",\
	 "normalized_ex_noise", "normalized_tot_noise"]
	for i in range (len(labels)):
		ws_spss.write(0, i, labels[i])
	labels = ["Genetic background", "normalized_in_noise_avg", "normalized_ex_noise_avg", \
	"normalized_tot_noise_avg", "nor_in_noise_err", "nor_ex_noise_err", "nor_tot_noise_err"]
	for i in range (len(labels)):
		ws_sum.write(0, i, labels[i])

	# read the standard background file to find the standard normalized 
	try:
		std_wb = xlrd.open_workbook(inputs[0], 'r')
	except XLRDError as e:
		print ('plot_normalized_five_backgrounds.py: Cannot open file"' + inputs[0] + '".')
		exit(1)
		
	# Find the average of standard mean
	std_ws = std_wb.sheet_by_name('Combined')
	file_len = std_ws.nrows
	std_noise = [[], [], []]
	std_mean = []
	for i in range (1, file_len): # slices
		row = list(std_ws.row(i))
		std_noise[intrinsic].append(row[in_input_index].value) # intrinsic noise
		std_noise[extrinsic].append(row[ex_input_index].value) # extrinsic noise
		std_noise[total].append(row[tot_input_index].value) # total noise, indices starting from 0
	for i in range(len(std_noise)):
		std_mean.append(numpy.mean(std_noise[i]))
	# END OF READING DATA OF STANDARD BACKGROUND FOR NORMALIZATION 
	
	#find and calculate data for each of the genetic backgrounds provided
	heights = [[], [], []] # in, ex, tot [in/ex/tot][background]
	errs = [[], [], []] # in, ex, tot
	row_index = 1
	for i in range (num_backgrounds):
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('plot_normalized_five_backgrounds.py: Cannot open file "'+inputs[i]+'".')
			exit(1)
			
		worksheet = workbook.sheet_by_name('Combined')	
		file_len =  worksheet.nrows
		normalized_noise = [[], [], []] # in, ex, tot
		for j in range(1, file_len): # slice						
			row = list(worksheet.row(j))
			# calculate and append necessary data into noises and errs
			(normalized_noise[intrinsic]).append(float(row[in_input_index].value) / float(std_mean[intrinsic]))
			(normalized_noise[extrinsic]).append(float(row[ex_input_index].value) / float(std_mean[extrinsic]))
			(normalized_noise[total]).append(float(row[tot_input_index].value) / float(std_mean[total])) # total noise / mean of the standard genetic backgrounds
			# write data
			ws_spss.write(row_index, 0, (i + 1))
			ws_spss.write(row_index, 1, row[in_input_index].value)
			ws_spss.write(row_index, 2, row[ex_input_index].value)
			ws_spss.write(row_index, 3, row[tot_input_index].value)
			ws_spss.write(row_index, 4, normalized_noise[intrinsic][j - 1])
			ws_spss.write(row_index, 5, normalized_noise[extrinsic][j - 1])
			ws_spss.write(row_index, 6, normalized_noise[total][j - 1])
			row_index += 1
			
		#append data of the bar_plots:
		for j in range(len(normalized_noise)):
			(heights[j]).append(numpy.mean(normalized_noise[j]))
			(errs[j]).append(2*numpy.std(normalized_noise[j])/math.sqrt(len(normalized_noise[j])))
		# write data
		line = [i + 1] + [heights[noise_index][i] for noise_index in range(len(heights))] + \
		[errs[noise_index][i] for noise_index in range(len(errs))]
		for k, item in enumerate(line):
			ws_sum.write(i + 1, k, item)
	return (heights,errs)

def plot_bar(num_geneticbackgrounds, ax, colors, names, heights, yerrs, ylabel):
	y_max = float('-inf')
	for i in range (num_geneticbackgrounds):
		this_height = heights[i] / heights[0] #normalize the average noise levels, based on the standard background, which is the first one in the input array
		y_max = max(y_max, this_height + 0.5 * yerrs[i])
		ax.bar(0.25 + i * (width + 0.3), this_height, width, yerr = yerrs[i], color=colors[i], ecolor='k', label = names[i]) 
		# draw normalized based on the first genetic background (wildtype or dmso)
	ax.set_yticks(numpy.arange(0, y_max+1, determineTickInterval(y_max,10)))
	ax.set_ylim(0,y_max*1.05)	
	ax.set_ylabel(ylabel, fontsize = 18)
	ax.tick_params(axis = 'x', which= 'both', bottom = 'off', top = 'off', labelbottom = 'off')
	ax.yaxis.set_tick_params(labelsize = 18)
	updateYTicklabels(ax)

def get_data_for_3noises_plot(heights3, errs3, heights2, errs2):
	# get rid of the standard genetic background noises
	heights = [[],[], []]
	errs = [[],[], []]
	for i in range(len(heights)):
		heights[i] = heights3[i][1:] + heights2[i][1:]
		errs[i] = errs3[i][1:] + errs2[i][1:]
	return heights, errs
	
def plot_3noises_3bg(heights, errs, directory, names, colors, num_bg):
	width = 0.5
	fig = plt.figure(figsize = (num_bg * 1.5 + width * (num_bg + 2), 6))
	ax = fig.add_subplot(111)
	for i in range(num_bg):#what genetic backgound
		xs = [(width * i + width * (j + 1) + width * num_bg * j) for j in range(len(heights))]
		ys = [heights[j][i] for j in range(len(heights))]
		yerrs = [errs[j][i] for j in range(len(errs))]
		ax.bar(xs, ys, width, color = colors[i], yerr = yerrs)
	xticks = [(width * (i + 1) + width * num_bg * i + width * num_bg * 0.5) for i in range (3)] # in, ex, tot
	ax.set_xlim(0, num_bg * 1.5 + width * (num_bg + 1))
	ax.set_xticks(xticks)
	ax.set_xticklabels(noise_name)
	# legend
	legend_shapes = []
	for i in range(len(names)):
		legend_shapes.append(mpatches.Patch(color=colors[i]))
	plt.legend(legend_shapes, names, loc = 9, bbox_to_anchor=(0.5,1.15), ncol=num_bg, numpoints=1, fontsize=14)
	# Adjust the subplot alignment
	fig.subplots_adjust(left=0.1, bottom=0.08 , right=0.95, top=0.85, wspace=None, hspace=0.3)
	plt.savefig(directory + "/normalized_all_noises.png", format = "png")
	
def compare_noise(inputs, names, colors, output_directory, num_geneticbackgrounds):
	for f in inputs:
		if not os.path.isfile(f):
			print ("plot_normalized_five_backgrounds.py: File '"+f+"' does not exist.")
			exit(1)

	# create plot
	f = plt.figure(figsize=(num_geneticbackgrounds * 1.25 + 2,6))
	gs = gridspec.GridSpec(1,2, width_ratios = [3,2])
	ax3 = f.add_subplot(gs[0])
	ax2 = f.add_subplot(gs[1])
	
	#create excel file
	write_workbook = xlwt.Workbook(encoding="ascii")
	wt_ws_sum = write_workbook.add_sheet("WTdeltaCdeltaD_summary")
	wt_ws_spss = write_workbook.add_sheet("WTdeltaCdeltaD_spss")
	dmso_ws_sum = write_workbook.add_sheet("dmsoDAPT_summary")
	dmso_ws_spss = write_workbook.add_sheet("dmsoDAPT_spss")
	
	# draw plot for wiltype, deltac, deltad
	(heights3, yerrs3) =  calculate_noise_data(3,inputs[:3], wt_ws_sum, wt_ws_spss) # get raw data of average noise and std
	tot_heights = [heights3[total][bg] for bg in range(len(heights3[0]))] # get total noise data
	tot_yerrs = [yerrs3[total][bg]  for bg in range(len(yerrs3[0]))] # get total errors
	plot_bar(3, ax3, colors[:3], names[:3], tot_heights, tot_yerrs, "Normalized Noise")
	
	# draw plot for dapt and dmso
	(heights2, yerrs2) = calculate_noise_data(2, inputs[3:5], dmso_ws_sum, dmso_ws_spss)
	tot_heights = [heights2[total][bg] for bg in range(len(heights2[0]))] # get total noise data
	tot_yerrs = [yerrs2[total][bg]  for bg in range(len(yerrs2[0]))] # get total errors
	plot_bar(2, ax2, colors[3:5],names[3:5], tot_heights, tot_yerrs, "Normalized Noise")

	
	# Figure legend
	legend_shapes = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(mpatches.Patch(color=colors[i]))
	plt.legend(legend_shapes, names, loc=0, bbox_to_anchor=(0.8,1.2), ncol=num_geneticbackgrounds, numpoints=1, fontsize=12)

	# Save figure
	f.subplots_adjust(left=0.1, bottom=0.1, right=.975, top=.80,  wspace=0.3, hspace=None)
	f.savefig(output_directory + "/normalized_noise_5bg.png", format = "png", dpi=300)
	'''
	# draw plots of all normalized noises for all genetic backgrounds, only plot the noise levels of deltac, deltad and dapt, normalized by its own standard background
	heights,errs = get_data_for_3noises_plot(heights3, yerrs3, heights2, yerrs2)
	plot_3noises_3bg(heights, errs, output_directory, names[1:3] + names[4:5], colors[1:3] + colors[4:5], 3)
	'''
	# save excel file
	write_workbook.save(output_directory + "/normalized_noise_five_backgrounds.xls")

def calculate_amplitude_data(num_bg, inputs, ws_sum, ws_spss):
	#write labels into the excel file
	labels = ["Genetic background", "average_normalized_amplitude", "std_error_amplitude"]
	for i in range (len(labels)):
		ws_sum.write(0, i, labels[i])
	labels = ["genetick background", "amplitude", "genetic background", "normalized_amp"]
	for i, label in enumerate(labels):
		ws_spss.write(0, i, label)
	# read the standard background file to find the standard normalized 
	try:
		std_wb = xlrd.open_workbook(inputs[0], 'r')
	except XLRDError as e:
		print ('plot_normalized_five_backgrounds.py: Cannot open file"' + inputs[0] + '".')
		exit(1)

	std_ws = std_wb.sheet_by_name('Her')
	row = list(std_ws.row(0)) # first row
	for j in range(1,len(row)):
		if row[j].value=='Average': # record mean amplitude
			std_avg_amp = std_ws.cell(1, j).value

	amp_her = [] # 1D: [background]
	std_her = [] # std error of the amplitude of this mutants, normalized, 1D: [background]
	row_index = 1
	for i in range(num_bg): # genetic background
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print ('compare_spatial_amplitude.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
			 						
		worksheet = workbook.sheet_by_name('Her') # her1 sheet	
		row = list(worksheet.row(0)) # first row
		nor_amp = []
		for j in range(1, len(row) - 1): # not looking at the first column and the last column, which contains meta data we do not need
			amp_raw = worksheet.cell(1, j).value
			amp_nor = worksheet.cell(1, j).value / std_avg_amp
			line = [i + 1, amp_raw, i + 1, amp_nor]
			#write data to spss
			for k, item in enumerate(line):
				ws_spss.write(row_index, k, item)
			row_index += 1
			nor_amp.append(worksheet.cell(1, j).value / std_avg_amp)
		mean = numpy.mean(nor_amp)
		ste = numpy.std(nor_amp) / math.sqrt(len(nor_amp))
		amp_her.append(mean)
		std_her.append(ste)
		sum_data = [i + 1, mean, ste]
		for j in range(len(sum_data)):
			ws_sum.write(i + 1, j, sum_data[j])
	return amp_her, std_her
		
def compare_amplitude(inputs, names, colors, output_directory, num_geneticbackgrounds):
	for f in inputs:
		if not os.path.isfile(f):
			print ("plot_normalized_five_backgrounds.py: File '"+f+"' does not exist.")
			exit(1)

	# create plot
	f = plt.figure(figsize=(num_geneticbackgrounds * 1.25 + 2,6))
	gs = gridspec.GridSpec(1,2, width_ratios = [3,2])
	ax3 = f.add_subplot(gs[0])
	ax2 = f.add_subplot(gs[1])

	#create excel file
	write_workbook = xlwt.Workbook(encoding="ascii")
	wt_ws_sum = write_workbook.add_sheet("WTdeltaCdeltaD_summary")
	wt_ws_spss = write_workbook.add_sheet("WTdeltaCdeltaD_spss")
	dmso_ws_sum = write_workbook.add_sheet("dmsoDAPT_summary")
	dmso_ws_spss = write_workbook.add_sheet("dmsoDAPT_spss")

	# draw plot for wiltype, deltac, deltad
	(heights3, yerrs3) =  calculate_amplitude_data(3,inputs[:3], wt_ws_sum, wt_ws_spss) # get raw data of average noise and std
	plot_bar(3, ax3, colors[:3], names[:3], heights3, yerrs3, "Total $\mathit{her}$ amplitude")
	
	(heights2, yerrs2) = calculate_amplitude_data(2, inputs[3:5], dmso_ws_sum, dmso_ws_spss)
	plot_bar(2, ax2, colors[3:5], names[3:5], heights2, yerrs2, "Total $\mathit{her}$ amplitude")
	
	# Figure legend
	legend_shapes = []
	for i in range(num_geneticbackgrounds):
		legend_shapes.append(mpatches.Patch(color=colors[i]))
	plt.legend(legend_shapes, names, loc=0, bbox_to_anchor=(0.8,1.2), ncol=num_geneticbackgrounds, numpoints=1, fontsize=12)

	# Save figure
	f.subplots_adjust(left=0.1, bottom=0.1, right=.975, top=.80,  wspace=0.3, hspace=None)
	f.savefig(output_directory + "/normalized_herAmp_5bg.png", format = "png", dpi=300)
	
	#save excel file
	write_workbook.save(output_directory + "/normalized_herAmp_five_backgrounds.xls")
	
def main():
	print ("Plotting normalized noise and amplitude for all five backgrounds....")
	# Check input
	if not shared.isInt(sys.argv[1]):
		print ('compare_normalized_five_backgrounds.py: Number of genetic backgrounds must be an integer.')
		exit(1)		
	num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine
	output_directory = sys.argv[2]
	if len(sys.argv)==3*num_geneticbackgrounds+3:
		inputs = sys.argv[3:num_geneticbackgrounds+3] # raw_noise.xls files
		names = sys.argv[num_geneticbackgrounds+3:2*num_geneticbackgrounds+3] # genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+3:] # colors for plotting different genetic backgrounds	
	else:
		usage()	
	
	assert num_geneticbackgrounds == 5
	compare_noise([directory + "/raw_noise.xls" for directory in inputs], names, colors, output_directory, num_geneticbackgrounds)
	compare_amplitude([directory + "/spatial_amplitude.xls" for directory in inputs], names, colors, output_directory, num_geneticbackgrounds)
	print ("Done.")
	
def usage():
	print ("compare_geneticbackgrounds.py: Invalid command-line arguments.")
	print ("Format: python compare_geneticbackgrounds.py <number of genetic backgrounds> <output directory> <path to input files for each genetic background> <genetic background names for labeling> <genetic background colors for plotting>")
	print ("Example: python plot_normalized_five_backgrounds.py 5 ../compare_output ../wildtypefulldataset/output/ ../deltacfulldataset/output/ ../deltadfulldataset/output/ ../DMSO/output/ ../DAPT/output/ Wildtype DeltaC DeltaD DMSO DAPT '#722AFF' g r '#FA5858' '#0000FF'")
	exit(1)

main()
