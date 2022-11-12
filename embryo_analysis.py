"""
Analyze experimental RNA expression data for a given embryo and write necessary Excel sheets for futher analyses
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
import matplotlib.pyplot as plt
import numpy, math
import xlrd, xlwt
from regions import Region
from xlrd import XLRDError

class Cell: # store the position and two expression levels for a cell
	def __init__(self, xpos, ypos, zpos, her1, her7, CB, YB):
		self.xpos = -xpos
		self.ypos = ypos
		self.zpos = zpos
		self.her1 = her1
		self.her7 = her7
		self.her1_bgN = her1 - CB	# Background normalization for her1 count
		self.her7_bgN = her7 - YB	# Background normalization for her7 count

def middle_splitting(cell_list): # split middle section into left(upper) and right(lower) sections
	ys = []
	for i in range(len(cell_list)):
		ys.append(cell_list[i].ypos)
	middle = (max(ys)+min(ys))/2
	upper = []
	lower = []
	for cell in cell_list:
		if cell.ypos >= middle:
			upper.append(cell)
		else: 
			lower.append(cell)
	return (upper, lower)

def analyze_slice(directory, wb, region): # extract necessary data from each slice, writes the data to slices.xls	
	# Set up the worksheet
	ws = wb.add_sheet("Region " + region.name)	
	labels = ["Slice #","# of cells","Her1 mean","Her1 variance","Her1 std","Her7 mean","Her7 variance","Her7 std"]

	for i in range(len(labels)):
		ws.write(0,i,label=labels[i])
		
	# Extract necessary data from each slice 
	for i in range(region.num_slices):
					
		curr_slice = region.slices[i]	
					 
		# Ignore slices with fewer than 3 cells
		if not curr_slice.identify_cells(region.cell_list):
			ws.write(i+1, 0, i+1)
			ws.write(i+1, 1, "Too few cells to analyze")
			continue
			
		mean_her1 = curr_slice.slice_mean_her1_bgN()
		mean_her7 = curr_slice.slice_mean_her7_bgN()
		var_her1 = curr_slice.slice_variance_her1_bgN()
		var_her7 = curr_slice.slice_variance_her7_bgN()
		std_her1 = numpy.sqrt(curr_slice.slice_variance_her1_bgN()) 
		std_her7 = numpy.sqrt(curr_slice.slice_variance_her7_bgN())	
				
		column_num = 0 
		line = [i+1,curr_slice.num_keep_cells,mean_her1,var_her1,std_her1,mean_her7,var_her7,std_her7]
		while column_num < len(line):
			ws.write(i+1, column_num, line[column_num])
			column_num+=1
		
		# Write individual cell expression levels
		for j in range(curr_slice.num_keep_cells):			
			ws.write(i+1, column_num, curr_slice.her1_bgNlevels[j])
			ws.write(i+1, column_num+1, curr_slice.her7_bgNlevels[j])
			column_num+=2

def write_slice_info(directory, wb, region):  # Write slices info to 'sliceInfo.xls' for heatmap plotting purpose

	ws_slices = wb.add_sheet("Slice " + region.name)
	labels = ["Slice #","bottom_left_xpos","bottom_right_xpos","top_left_pos","top_right_pos","# of cells","Her1 mean","Her1 variance","Her1 std","Her7 mean","Her7 variance","Her7 std"]
	
	for i in range(len(labels)):
		ws_slices.write(3,i,label=labels[i])
		
	ws_slices.write(0,1,label="top")
	ws_slices.write(0,2,label="bottom")

	for i in range(region.num_slices):					
		curr_slice = region.slices[i]
		
		if i == 0:
			ws_slices.write(1,1,curr_slice.top)
			ws_slices.write(1,2,curr_slice.bottom)
				
		ws_slices.write(i+4,0,i+1)
		ws_slices.write(i+4,1,curr_slice.bottom_left_xpos)
		ws_slices.write(i+4,2,curr_slice.bottom_right_xpos)
		ws_slices.write(i+4,3,curr_slice.top_left_xpos)
		ws_slices.write(i+4,4,curr_slice.top_right_xpos)
		ws_slices.write(i+4,5,curr_slice.num_cells)

		# Ignore slices with fewer than 3 cells
		if curr_slice.num_cells <= 2:
			ws_slices.write(i+4, 6, "Too few cells to analyze")
			continue
		
		ws_slices.write(i+4,6,curr_slice.slice_mean_her1())
		ws_slices.write(i+4,7,curr_slice.slice_variance_her1())
		ws_slices.write(i+4,8,numpy.sqrt(curr_slice.slice_variance_her1()) )
		ws_slices.write(i+4,9,curr_slice.slice_mean_her7())
		ws_slices.write(i+4,10,curr_slice.slice_variance_her7())
		ws_slices.write(i+4,11,numpy.sqrt(curr_slice.slice_variance_her7()))
		
		column_num = 12
		# Write individual cell expression levels
		for j in range(curr_slice.num_cells):
			ws_slices.write(i+4, column_num, curr_slice.her1_levels[j])
			ws_slices.write(i+4, column_num+1, curr_slice.her7_levels[j])
			column_num+=2
	
def plother1her7(cell_lists, directory):
	her1 = []
	her7 = []
	her = []
	for i in range(len(cell_lists)):
		for j in range(len(cell_lists[i])):
			her1.append(cell_lists[i][j].her1)
			her7.append(cell_lists[i][j].her7)
			her.append(cell_lists[i][j].her1 + cell_lists[i][j].her7)

	# Plot her1 & her7 expression distribution
	plt.subplot(211)
	#plt.hist(her7, 50, normed=1, facecolor='lightblue', alpha=1.0, label='her7')
	#plt.hist(her1, 50, normed=1, facecolor='orange', alpha=0.5, label='her1')
	plt.hist(her7, 50, facecolor='lightblue', alpha=1.0, label='her7')
	plt.hist(her1, 50, facecolor='orange', alpha=0.5, label='her1')
	plt.legend()
	plt.tick_params(direction='in')
	plt.ylabel('Frequency')
	#plt.xlabel('mRNA expression level')
	#plt.savefig(directory + "/her1her7hist.png", format = "png")

	# Plot total her expression distribution
	plt.subplot(212)
	#plt.hist(her, 50, normed=1, facecolor='orange', alpha=1.0, label='Total her')	
	plt.hist(her, 50, facecolor='orange', alpha=1.0, label='Total her')
	plt.legend()
	plt.tick_params(direction='in')
	plt.ylabel('Frequency')
	plt.xlabel('mRNA expression level')
	plt.savefig(directory + "/totalherhist.png", format = "png", dpi=300)
		
def main():
	args = sys.argv[1:]
	num_args = len(args)
	req_args = [False] * 8
	lr_shift = 0
	ly_shift = 0
	middle = 0
	wholePSM = False
	if num_args >= 16:
		for arg in range(0, num_args - 1, 2):
			option = args[arg]
			value = args[arg + 1]
			
			if option == '-i' or option == '--input-file':
				filename = value				
				if not os.path.isfile(filename):
					print("embryo_analysis.py: File "+filename+" does not exist.")
					exit(1)
				req_args[0] = True
			elif option == '-d' or option == '--output-directory': # output directory
				directory = value
				req_args[1] = True
			elif (option == '-a' or option == '--initial-angle') and shared.isFloat(value):	# y-intercept from angle measurement equation
				L_angle = 180 - float(value)
				R_angle = 180 + float(value)
				req_args[2] = True
			elif (option == '-dA' or option == '--delta-angle') and shared.isFloat(value): # Slope from angle measurement equation, use 0.0 for fix angle
				L_delta_angle = -float(value)	# angle decrease after every step in left PSM
				R_delta_angle = float(value)	# angle increase after every step in right PSM
				req_args[3] = True
			elif (option == '-n' or option == '--num-sec') and shared.isInt(value):
				num_sec = int(value)
				req_args[4] = True				
			elif (option == '-m1' or option == '--background-noise-mean-her1') and shared.isFloat(value):
				CB = float(value)
				req_args[5] = True
			elif (option == '-m7' or option == '--background-noise-mean-her7') and shared.isFloat(value):
				YB = float(value)
				req_args[6] = True
			elif (option == '-f' or option == '--input-format') and shared.isInt(value): # whether the data already split in half and shifted
				in_format = True
				if int(value)==0: # 1 if the data is already split and shifted, 0 if the data needs to be split and shifted
					in_format = False
				req_args[7] = True
				
			# (Optional) If the data needs to be split not exactly in half, specify how much the threshold for left and right should be shifted
			elif (option == '-s' or option == '--shift') and shared.isFloat(value): 
				lr_shift = float(value)
			# (Optional) If the image is from lateral view without left and right PSM info
			elif (option == '-w' or option == '--wholePSM') and shared.isInt(value): 
				if int(value)==1: # 1 if image is from lateral view
					wholePSM = True
			# (Optional) Some mutants require the left section to be shifted along the y-axis. This specifies how much shift is needed.
			elif (option == '-ly' or option == '--left-yaxis-shift') and shared.isFloat(value):
				ly_shift = float(value)
			# (Optional) Number of middle sections if they exist (wildtype only)
			elif (option == '-m' or option == '--middle-section') and shared.isInt(value):
				middle = int(value)
			# (Optional) Value for left angle - use only if left and right initiate angle is different
			elif (option == '-l' or option == '--l-angle') and shared.isFloat(value):
				L_angle = float(value)
			# (Optional) Value for right angle - use only if left and right initiate angle is different	
			elif (option == '-r' or option == '--r-angle') and shared.isFloat(value):
				R_angle = float(value)
			elif option == '-h' or option == '--help':
				usage()
			else:
				usage()
		for arg in req_args:
			if not arg:
				usage()
	else:
		usage()

	
	shared.ensureDir(directory)

	# Open the input file
	workbook = xlrd.open_workbook(filename,'r')
	worksheets = workbook.sheet_names()

	# Iterate over each worksheet in a workbook
	for worksheet_name in worksheets:
        
		worksheet = workbook.sheet_by_name(worksheet_name)
		file_len = worksheet.nrows
		
		cell_lists = []
		left_cell_lists = []
		right_cell_lists = []
		middle_cell_lists = []
		section_xmin = [0] * num_sec
		section_xmax = [0] * num_sec
		ymin = 0
		ymax = 0

		if wholePSM:	# If image is lateral view, without left/right partition
			print(filename)
			cell_lists.append([])
			for j in range(1,file_len):
				row = list(worksheet.row(j))
				for i in range(num_sec):
					if row[i*6].value != '' and isinstance(row[i*6+1].value,float):
						curr_cell = Cell(row[i*6+1].value, row[i*6+2].value+ly_shift, row[i*6+3].value, row[i*6+4].value, row[i*6+5].value, CB, YB)	
						cell_lists[0].append(curr_cell)
			
			regions = [Region(num_sec, cell_lists, "L", L_angle, L_delta_angle)]

			workbook = xlwt.Workbook(encoding="ascii")
			for i in range(len(regions)):
				worksheet = workbook.add_sheet("Region " + regions[i].name)	
				
				# Slice boundary data
				labels = ["Top left xpos","Top right xpos","Top ypos","Bottom left xpos","Bottom right xpos","Bottom ypos","Slice width","# of slices"]
				line = [regions[i].slices[0].top_left_xpos, regions[i].slices[-1].top_right_xpos, regions[i].top.ypos,
					regions[i].slices[0].bottom_left_xpos, regions[i].slices[-1].bottom_right_xpos, regions[i].bottom.ypos,
					regions[i].slice_width, regions[i].num_slices]		
				for j in range(len(labels)):
					worksheet.write(0, j, labels[j])
					worksheet.write(1, j, line[j])
							
				labels = ["Cell xpos","Cell ypos","Her1 level","Her7 level"]
				for j in range(len(labels)):
					worksheet.write(3,j,labels[j])
				for j in range(len(regions[i].cell_list)):
					worksheet.write(j+4, 0, regions[i].cell_list[j].xpos)
					worksheet.write(j+4, 1, regions[i].cell_list[j].ypos)
					worksheet.write(j+4, 2, regions[i].cell_list[j].her1)
					worksheet.write(j+4, 3, regions[i].cell_list[j].her7)
											
			workbook.save(directory + "/cells.xls")							
			workbook = xlwt.Workbook(encoding="ascii")		
			analyze_slice(directory, workbook, regions[0])
			workbook.save(directory + "/slices.xls")	# Write background normalized her info for every slices, for downstream analysis
			workbook = xlwt.Workbook(encoding="ascii")
			write_slice_info(directory, workbook, regions[0])
			workbook.save(directory + "/sliceInfo.xls")	# Write raw her count for every slice, for heatmap plotting and visualization purpose
		
		else:		# If image is in superior view, with left/right partition
				
			# Initialize the lists holding each section's list of cells
			for i in range(num_sec):
				cell_lists.append([])	
				left_cell_lists.append([])
				right_cell_lists.append([])
			
			# Put the data coming from the files to the matrix
			for j in range(1, file_len):
				row = list(worksheet.row(j))
				
				if in_format: # given data is already split in half and shifted
					for i in range(num_sec * 2):
						if row[i*6].value != '' and isinstance(row[i*6+1].value,float):
							cur_i = int(i/2)
							if i%2 == 0:
								curr_cell = Cell(row[i*6+1].value, row[i*6+2].value+ly_shift, row[i*6+3].value, row[i*6+4].value, row[i*6+5].value, CB, YB)	
								left_cell_lists[cur_i].append(curr_cell)
							else:
								curr_cell = Cell(row[i*6+1].value, row[i*6+2].value, row[i*6+3].value, row[i*6+4].value, row[i*6+5].value, CB, YB)	
								right_cell_lists[cur_i].append(curr_cell)
							cell_lists[cur_i].append(curr_cell)
					
					if middle!=0 and row[num_sec*2*6] != '' and isinstance(row[num_sec*2*6+1].value,float):
						curr_cell = Cell(row[num_sec*12+1].value, row[num_sec*12+2].value+ly_shift, row[num_sec*12+3].value, row[num_sec*12+4].value, row[num_sec*12+5].value, CB, YB)	
						middle_cell_lists.append(curr_cell)
				
				else: # given data will be split in half and shifted later			
					for i in range(num_sec):
						if row[i*6].value != '':
							curr_cell = Cell(row[i*6+1].value, row[i*6+2].value, row[i*6+3].value, row[i*6+4].value, row[i*6+5].value, CB, YB)
							cell_lists[i].append(curr_cell)
							if section_xmin[i] == 0 or section_xmin[i] > row[i*6+1].value:
								section_xmin[i] = row[i*6+1].value
							if section_xmax[i] == 0 or section_xmax[i] < row[i*6+1].value:
								section_xmax[i] = row[i*6+1].value
							if ymin == 0 or ymin > row[i*6+2].value:
								ymin = row[i*6+2].value
							if ymax == 0 or ymax < row[i*6+2].value:
								ymax = row[i*6+2].value

			if in_format and len(middle_cell_lists)>0:
				upper, lower = middle_splitting(middle_cell_lists)
				left_cell_lists.append(upper)
				right_cell_lists.append(lower)
								
			if not in_format:
				# Shift sections to correct positions
				shift = 0
				for i in range(num_sec-1,0,-1):
					shift += section_xmax[i]-section_xmin[i-1] + 0.01
					for j in range(len(cell_lists[i-1])):
						cell_lists[i-1][j].xpos += shift

				# Split sections into left and right
				split_threshold = (ymax+ymin)/2 + lr_shift # shift the threshold up or down (left or right)
				for i in range(num_sec):
					for j in range(len(cell_lists[i])):
						if cell_lists[i][j].ypos > split_threshold:
							left_cell_lists[i].append(cell_lists[i][j])
						else:
							right_cell_lists[i].append(cell_lists[i][j])

			regions = [Region(num_sec, left_cell_lists, "L", L_angle, L_delta_angle), Region(num_sec, right_cell_lists, "R", R_angle, R_delta_angle)]
			
			workbook = xlwt.Workbook(encoding="ascii")
			for i in range(len(regions)):
				worksheet = workbook.add_sheet("Region " + regions[i].name)	
				
				# Slice boundary data
				labels = ["Top left xpos","Top right xpos","Top ypos","Bottom left xpos","Bottom right xpos","Bottom ypos","Slice width","# of slices"]
				line = [regions[i].slices[0].top_left_xpos, regions[i].slices[-1].top_right_xpos, regions[i].top.ypos,
					regions[i].slices[0].bottom_left_xpos, regions[i].slices[-1].bottom_right_xpos, regions[i].bottom.ypos,
					regions[i].slice_width, regions[i].num_slices]		
				for j in range(len(labels)):
					worksheet.write(0, j, labels[j])
					worksheet.write(1, j, line[j])
							
				labels = ["Cell xpos","Cell ypos","Her1 level","Her7 level"]
				for j in range(len(labels)):
					worksheet.write(3,j,labels[j])
				for j in range(len(regions[i].cell_list)):
					worksheet.write(j+4, 0, regions[i].cell_list[j].xpos)
					worksheet.write(j+4, 1, regions[i].cell_list[j].ypos)
					worksheet.write(j+4, 2, regions[i].cell_list[j].her1)
					worksheet.write(j+4, 3, regions[i].cell_list[j].her7)
											
			workbook.save(directory + "/cells.xls")							
			workbook = xlwt.Workbook(encoding="ascii")		
			analyze_slice(directory, workbook, regions[0])
			analyze_slice(directory, workbook, regions[1])
			workbook.save(directory + "/slices.xls")	# Write background normalized her info for every slices, for downstream analysis
			workbook = xlwt.Workbook(encoding="ascii")
			write_slice_info(directory, workbook, regions[0])
			write_slice_info(directory, workbook, regions[1])
			workbook.save(directory + "/sliceInfo.xls")	# Write raw her count for every slice, for heatmap plotting and visualization purpose
			
		# make histogram to view distribution of her1/her7 expression level
		plother1her7(cell_lists,directory)
		break # force the for loop to run only once (read the first worksheet only)

def usage():
	print("embryo_analysis.py: Invalid command-line arguments.")
	print("Format: python embryo_analysis.py -i <input Excel file> -d <output directory> -a <initial angle from posterior> -dA <angle change rate> -n <number of sections> -m1 <background-noise-mean-her1> -m7 <background-noise-mean-her7> -f <0 or 1 to specify input format> -s <optional:half threshold shift> -l <optional:angle for left PSM> -r <optional:angle for right PSM>")
	print("Example: python embryo_analysis.py -i wildtypefulldataset/WT1.xlsx -d wildtypefulldataset/embryo1 -a 44.23 -dA 0.039 -n 6 -m1 0.019 -m2 0.076 -f 0 -s -20")
	exit(1)

main()
