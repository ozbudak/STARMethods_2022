'''
Combine spatial expression data and calculate spatial amplitude of oscillations in gene expression
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
'''
import sys, shared, os
import numpy, math
import matplotlib.pyplot as plt
import xlrd, xlwt
from xlrd import XLRDError
from scipy.ndimage import gaussian_filter1d # Gaussian filtering package
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

colors = ['b','r','g'] # colors for her1, her7, and her
width = 0.4 # bar graph width

default_num_slice_per_section = 5
num_slice_per_section = 5
default_percent_top_bottom_amp = 10
percent_top_bottom_amp = 10

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.1,0.5,1,2,5,10,20,30,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def smoothen(x): # apply Gaussian smoothing to a given vector
	window = 3 # window size
	sigma = 10 # sigma
	for idx in range(0, len(x) - window + 1):
		smooth_slice = gaussian_filter1d(x[idx:idx+window], sigma)
		for idxx in range(len(smooth_slice)):
			x[idx + idxx] = smooth_slice[idxx]	
	return x # return smoothened x vector

def findAmplitudeSpace(left_slices, right_slices, num_slices): # measure amplitudes
	amp = []
	stderr = []
	for i in range(int(num_slices/num_slice_per_section)): # every five spatial locations (cell positions)
		s = []
		for j in range(i*num_slice_per_section,(i+1)*num_slice_per_section):
			s += left_slices[j]
			if len(right_slices)>0:
				s += right_slices[j]		
		s = [x for x in s if not math.isnan(x)]
		s = sorted(s) # rank all slices wihtin the five spatial locations
		#indices for top and bottom indices
		bottom_index = int(len(s) * float(percent_top_bottom_amp) / float(100))
		top_index = int(len(s) * float(100 - percent_top_bottom_amp) / float(100))
		bottom = s[:bottom_index] # bottom percent_top_bottom_amp%
		top = s[top_index:] # top percent_top_bottom_amp%
		
		# Amplitude = difference between the averages of top and bottom
		amp.append(numpy.mean(top)-numpy.mean(bottom))
				
		# Standard deviation = sqrt(std(top)^2+std(bottom)^2)
		# Standard error = standard deviation / sqrt(|top|+|bottom|)
		stderr.append(math.sqrt(numpy.std(top)**2 + numpy.std(bottom)**2) / math.sqrt(len(top+bottom)))					
	
	return (amp,stderr)	

def interpolate(x): # interpolate empty data points
	for i in range(1, len(x)-1): # looking for one empty data point
		if numpy.isnan(x[i]) and not numpy.isnan(x[i-1]) and not numpy.isnan(x[i+1]):
			x[i] = (x[i-1] + x[i+1])/2 # average of two neighbors
	for i in range(1, len(x)-2): # looking for two consecutive empty data points
		if numpy.isnan(x[i]) and numpy.isnan(x[i+1]) and not numpy.isnan(x[i-1]) and not numpy.isnan(x[i+2]):
			x[i] = (x[i-1] + x[i+2])/2 # average of two outer neighbors 
			x[i+1] = (x[i-1] + x[i+2])/2
	return x

def determineStartEnd(x): # determine where data starts/ends to be valid (mutant may start/end with slices with too few cells)
	i = 0
	while i < (len(x)) and numpy.isnan(x[i]):
		i+=1	
	j = len(x)-1
	while j > 0 and numpy.isnan(x[j]):
		j-=1	
	return i, j+1 # start, end

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()] # intergers
	ax.set_xticklabels(xlabels)
	updateYTicklabels(ax)
	
def updateYTicklabels(ax):
	ylabels = [format(label, r',.0f') for label in ax.get_yticks()] # intergers
	ax.set_yticklabels(ylabels)

def plotAmplitudeBar(ax, amplitude_space_her1, amplitude_space_her7): # plot bar graph showing amplitude in her1 and her7
	ax.bar(0.2,numpy.mean(amplitude_space_her1),width,yerr=2*numpy.std(amplitude_space_her1)/math.sqrt(len(amplitude_space_her1)),color='b',ecolor='k',capsize=4, linewidth=0.5)
	ax.bar(1.2,numpy.mean(amplitude_space_her7),width,yerr=2*numpy.std(amplitude_space_her7)/math.sqrt(len(amplitude_space_her7)),color='r',ecolor='k',capsize=4, linewidth=0.5)	
	
	ax.tick_params(direction='in')
	ax.set_xticks([0+width/2, 1+width/2])
	ax.set_xticklabels((r'\textit{her1}',r'\textit{her7}'))
	ax.tick_params(axis='x', pad=5, direction='in')
	ax.set_xlim(-0.3,1.6)	
	
	start, end = ax.get_ylim()
	ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
	ax.set_ylim(0,end)
	ax.set_ylabel(r'Spatial amplitude')
	updateYTicklabels(ax)

def plotAmplitudeSpace(axes, amplitude_space, amplitude_space_stderr): # plot amplitudes across space		
	x_loc = []
	for i in range(len(amplitude_space[0])):
		x_loc.append(i*5+2.5)
		
	for i in range(len(axes)):
		axes[i].plot(x_loc,amplitude_space[i], c=colors[i])
		for j in range(len(x_loc)):			
			axes[i].errorbar(x_loc[j], amplitude_space[i][j], yerr=2*amplitude_space_stderr[i][j], ls='none', c=colors[i], capsize=2)

	for ax in axes:
		start, end = ax.get_xlim()
		ax.set_xticks(numpy.arange(0, end+1, 10)) # x-axis tick interval set to 10		
		start, end = ax.get_ylim()
		ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
		ax.set_ylim(0,end)
		ax.tick_params(direction='in')
		updateTicklabels(ax)
	
	axes[2].set_xlabel('Cell position (posterior - anterior)')
	axes[0].set_ylabel(r'\textit{her1} amplitude')
	axes[1].set_ylabel(r'\textit{her7} amplitude')
	axes[2].set_ylabel(r'\textit{her} amplitude')

def plotCombinedScatter(axes, slices, num_slices): # plot combined expression data across space as a scatter plot
	for i in range(len(axes)): # axis
		for j in range(num_slices): # cell position
			axes[i].scatter([j]*len(slices[i][j]), slices[i][j], c=colors[i], s = 22, edgecolors='none')			
	axes[0].set_ylabel(r'\textit{her1} mRNA')
	axes[1].set_ylabel(r'\textit{her7} mRNA')
	axes[2].set_ylabel(r'\textit{her} mRNA')
	axes[2].set_xlabel('Cell position (posterior - anterior)')	
	for ax in axes:
		ax.set_xlim(-2,num_slices+1)
		start, end = ax.get_xlim()
		ax.set_xticks(numpy.arange(0, end+1, 10)) # x-axis tick interval set to 10		
		start, end = ax.get_ylim()
		ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
		ax.set_ylim(0,end)
		ax.tick_params(direction='in')
		updateTicklabels(ax)
		
def plotCombinedLine(axes, slices, num_slices): # plot combined expression data across space as a line graph	
	for i in range(len(axes)):
		for j in range(num_slices):
			slices[i][j] = [x for x in slices[i][j] if not math.isnan(x)]
			axes[i].scatter(j, numpy.mean(slices[i][j]), c=colors[i], s = 22, edgecolors='none')
			axes[i].errorbar(j, numpy.mean(slices[i][j]), yerr=2*numpy.std(slices[i][j])/math.sqrt(len(slices[i][j])), ls='none', c=colors[i],capsize=2)
	axes[2].set_xlabel('Cell position (posterior - anterior)')
	for ax in axes:
		ax.set_xlim(-2,num_slices+1)
		start, end = ax.get_xlim()
		ax.set_xticks(numpy.arange(0, end+1, 10)) # x-axis tick interval set to 10		
		start, end = ax.get_ylim()
		ax.set_yticks(numpy.arange(0, end+1, determineTickInterval(end,6)))
		ax.set_ylim(0,end)
		ax.tick_params(direction='in')
		updateTicklabels(ax)

def writeAmplitudeData(ws, amplitude, amplitude_stderr): # write amplitude data
	ws.write(0,0,'Group # (position)')
	ws.write(1,0,'Amplitude')
	ws.write(2,0,'Std error')	
	for i in range(len(amplitude)):
		ws.write(0,i+1,str(i+1)+' ('+str(i*5)+'-'+str(i*5+4)+')') # five cell positions
		ws.write(1,i+1,amplitude[i])
		ws.write(2,i+1,amplitude_stderr[i])	
	ws.write(0,len(amplitude)+1,'Average')
	ws.write(1,len(amplitude)+1,numpy.mean(amplitude))
	ws.write(2,len(amplitude)+1,numpy.std(amplitude)/math.sqrt(len(amplitude)))

def writeCombinedData(ws, left_slices, right_slices, num_embryos, num_slices): # write combined data
	labels = ['Cell position','Mean','Stdev','Min','Max']
	for i in range(len(labels)):
		ws.write(0,i,labels[i])
	for i in range(num_embryos):
		ws.write(0,i+len(labels),'Embryo '+str(i+1)+' L')
		ws.write(0,i+num_embryos+len(labels),'Embryo '+str(i+1)+' R')	
	for i in range(num_slices):		
		tmp_left = [x for x in left_slices[i] if not math.isnan(x)] # temporarily get rid of nan values to compute mean and standard deviation
		tmp_right = []
		if len(right_slices)>0:
			tmp_right = [x for x in right_slices[i] if not math.isnan(x)]
		combined = tmp_left+tmp_right
		if len(combined)>0:
			line = [i,numpy.mean(combined),numpy.std(combined),min(combined),max(combined)]		
			for j in range(len(line)):
				ws.write(i+1,j,line[j])
		for k in range(len(left_slices[i])): # all left slices
			if numpy.isnan(left_slices[i][k]): # blank if data invalid
				left_slices[i][k]=''
			ws.write(i+1,k+5,left_slices[i][k])
			
		if len(right_slices)==0:	# skip for lateral view PSM comparison
			continue
			
		for k in range(len(right_slices[i])): # all right slices
			if numpy.isnan(right_slices[i][k]): # blank if data invalid
				right_slices[i][k]=''
			ws.write(i+1,k+5+len(left_slices[i]),right_slices[i][k]) # right slice data written after all left slice data
	
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print('plot_spatial_amplitude.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print('plot_spatial_amplitude.py: Number of embryos must be larger than zero.')
		exit(1)		
	
	num_embryos = int(sys.argv[1])
	    	        
	if len(sys.argv)==num_embryos+3:
		inputs = sys.argv[2:num_embryos+2]
		directory = sys.argv[num_embryos+2]	
	else:
		usage()
		        
	shared.ensureDir(directory) 
	
	left_her1 = []
	left_her7 = []
	left_her = []
	right_her1 = []
	right_her7 = []
	right_her = []
	num_slices = []        
                
	for i in range(len(inputs)): # embryo
		# Open embryo data
		if not os.path.isfile(inputs[i]):
			print('plot_spatial_amplitude.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print('plot_spatial_amplitude.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		worksheets = workbook.sheets()	
			
		for j in range(len(worksheets)): # region
			worksheet = worksheets[j]
			file_len =  worksheet.nrows
			num_slices.append(file_len-1)
						
			for k in range(1,file_len): # slice (cell position)
				row = list(worksheet.row(k)) # read data
				
				if j==0 and len(left_her1)<k: # two-dimensional arrays
					left_her1.append([])					
					left_her7.append([])	
					left_her.append([])			
				elif j==1 and len(right_her1)<k:
					right_her1.append([])
					right_her7.append([])
					right_her.append([])					
							
				if isinstance(row[1].value,float): # enough cells in this slice to analyze
					num_cells = int(row[1].value)
					her1 = [] # expression level from each cell
					her7 = []
					
					for l in range(num_cells): 
						h1 = row[8+2*l].value # cell's her1 expression level
						h7 = row[8+2*l+1].value # cell's her7 expression level
						
						# Set negative expression levels to zero
						her1.append(h1 if h1 > 0 else 0) 
						her7.append(h7 if h7 > 0 else 0)					
													
					her1_mean = numpy.mean(her1)
					her7_mean = numpy.mean(her7)
					her_mean = her1_mean+her7_mean	
											
					if j==0: # left	
						left_her1[k-1].append(her1_mean)
						left_her7[k-1].append(her7_mean)
						left_her[k-1].append(her_mean)							
					else: # right					
						right_her1[k-1].append(her1_mean)
						right_her7[k-1].append(her7_mean)
						right_her[k-1].append(her_mean)	
					
				else: # too few cells to analyze				
					if j==0:
						left_her1[k-1].append(numpy.nan) # nan as placeholders
						left_her7[k-1].append(numpy.nan)
						left_her[k-1].append(numpy.nan)
					else:
						right_her1[k-1].append(numpy.nan)
						right_her7[k-1].append(numpy.nan)
						right_her[k-1].append(numpy.nan)		
	
	
	# Determine number of slices found in at least 80% of embryos for analysis
	nSlices = sorted(num_slices)[int(num_embryos*0.2)]
	
	### Measure and plot amplitude for every five cell positions ###
	amplitude_space_her1,amplitude_space_stderr_her1 = findAmplitudeSpace(left_her1, right_her1, nSlices)
	amplitude_space_her7,amplitude_space_stderr_her7 = findAmplitudeSpace(left_her7, right_her7, nSlices)
	amplitude_space_her,amplitude_space_stderr_her = findAmplitudeSpace(left_her, right_her, nSlices)
	
	# Bar graph showing her1 and her7 average amplitudes (does not show how amplitude changes in space)
	fig = plt.figure(figsize=(3,4),dpi=300)	
	ax = fig.add_subplot(111)
	plotAmplitudeBar(ax, amplitude_space_her1, amplitude_space_her7)
	fig.subplots_adjust(left=0.2, bottom=0.075, right=0.95, top=.95,  wspace=None, hspace=0.2)	
	fig.savefig(directory + '/average_spatial_amplitude_bar.png', format = 'png', dpi=300)
	
	# Line graph showing amplitude across space
	fig = plt.figure(figsize=(5,9),dpi=300)
	axes = [fig.add_subplot(311),fig.add_subplot(312),fig.add_subplot(313)]		
	plotAmplitudeSpace(axes, [amplitude_space_her1, amplitude_space_her7, amplitude_space_her], [amplitude_space_stderr_her1, amplitude_space_stderr_her7, amplitude_space_stderr_her])
	fig.subplots_adjust(left=0.15, bottom=0.075, right=0.95, top=.975,  wspace=None, hspace=0.2)	
	fig.savefig(directory + '/spatial_amplitude_line.png', format = 'png', dpi=300)
	
	# Write amplitude data
	workbook = xlwt.Workbook(encoding='ascii') # output Excel file
	writeAmplitudeData(workbook.add_sheet('Her1'), amplitude_space_her1, amplitude_space_stderr_her1)
	writeAmplitudeData(workbook.add_sheet('Her7'), amplitude_space_her7, amplitude_space_stderr_her7)
	writeAmplitudeData(workbook.add_sheet('Her'), amplitude_space_her, amplitude_space_stderr_her)	
	workbook.save(directory + '/spatial_amplitude.xls')
	
	
	
	# Combine spatial expression data from left and right regions
	slices_her1 = left_her1+right_her1
	slices_her7 = left_her7+right_her7
	slices_her = left_her+right_her
	
	### Plot combined spatial expression data ###	
	fig = plt.figure(figsize=(9,9))	
	axes = [fig.add_subplot(321), fig.add_subplot(323), fig.add_subplot(325)]
	plotCombinedScatter(axes, [slices_her1, slices_her7, slices_her], nSlices) # scatter plot on the left plotting all data points	
			
	axes = [fig.add_subplot(322),fig.add_subplot(324),fig.add_subplot(326)]
	plotCombinedLine(axes, [slices_her1, slices_her7, slices_her], nSlices) # line graph on the right showing average and two standard errors	
	fig.subplots_adjust(left=0.1, bottom=None, right=0.975, top=.95,  wspace=None, hspace=0.2)
	fig.savefig(directory + '/combined_spatial_expression.png', format = 'png')
	
	### Create plots as one of the 
	# Write combined spatial expression data to Excel
	workbook = xlwt.Workbook(encoding='ascii')
	writeCombinedData(workbook.add_sheet('Her1'),left_her1,right_her1,len(inputs),nSlices)
	writeCombinedData(workbook.add_sheet('Her7'),left_her7,right_her7,len(inputs),nSlices)
	writeCombinedData(workbook.add_sheet('Her'),left_her,right_her,len(inputs),nSlices)	
	workbook.save(directory + '/combined_spatial_expression.xls')
	
def usage():
	print('plot_spatial_amplitude.py: Invalid command-line arguments')
	print("Format: python plot_spatial_amplitude.py <number of embryos> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print('Example: python plot_spatial_amplitude.py 20 ../wildtypefulldataset/output/embryo1/slices.xls ../wildtypefulldataset/output/embryo2/slices.xls ... ../wildtypefulldataset/output/embryo20/slices.xls ../wildtypefulldataset')
	exit(1)

main()
