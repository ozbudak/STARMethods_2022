"""
Create two-dimensional heatmaps consisting of slice boundaries and individual cells
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
import math
import matplotlib.pyplot as plt
import sys, shared, os
import xlrd
from xlrd import XLRDError

lw = 0.2 # line width for slice boundaries
cs = 40 # cell size for plotting 
colorred = ["#f2f3f4","#fdedec","#f1948a","#ec7063","#e74c3c","#cb4335","#b03a2e","#943126","#641e16","#2e4053"]
colorblue = ["#f2f3f4","#ebf5fb","#85c1e9","#5dade2","#3498db","#2e86c1","#2874a6","#21618c","#1b4f72","#2e4053"]

def main():
	args = sys.argv[1:]
	num_args = len(args)
	req_args = [False]*3
	
	if num_args >= 4:
		for arg in range(0, num_args - 1, 2):
			option = args[arg]
			value = args[arg + 1]
			if option == '-i' or option == '--input-file': # cells.xls
				filename = value
				if not os.path.isfile(filename):
					print("create_heatmap.py: File "+filename+" does not exist.")
					exit(1)
				req_args[0] = True
			elif option == '-s' or option == '--slice-Info': # sliceInfo.xls
				sliceinfo_file = value
				if not os.path.isfile(sliceinfo_file):
					print("create_heatmap.py: File "+sliceinfo_file+" does not exist.")
					exit(1)
				req_args[1] = True
			elif option == '-d' or option == '--output-directory':
				directory = value
				req_args[2] = True			
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
	
	try:
		workbook = xlrd.open_workbook(filename,'r')
	except XLRDError as e:
		print("create_heatmap.py: Cannot open "+filename+".")
		exit(1)	
	worksheets = workbook.sheets()

	try:
		workbook_slice = xlrd.open_workbook(sliceinfo_file,'r')
	except XLRDError as e:
		print("create_heatmap.py: Cannot open "+sliceinfo_file+".")
		exit(1)
		
	# Set up two figures (her1 and her7) --plotted simultaneously
	fig_her1 = plt.figure(figsize=(8,4))	
	ax_her1 = fig_her1.add_subplot(111)
	fig_her7 = plt.figure(figsize=(8,4))
	ax_her7 = fig_her7.add_subplot(111)
	
	top = 0 # y coordinate of the top end of the embryo
	bottom = 0 # y coordinate of the bottom end of the embryo
	left = float("inf") # x coordinate of the left end of the embryo
	right = float("-inf") # x coordinate of the right end of the embryo

	# Plot each region (left and right)
	for i in range(len(worksheets)):								
		# Read slices info
		if i == 0:
			worksheet_slice = workbook_slice.sheet_by_name("Slice L")
		else:
			worksheet_slice = workbook_slice.sheet_by_name("Slice R")
			
		file_len_slice = worksheet_slice.nrows
		
		row = list(worksheet_slice.row(1))
		top_ypos = row[1].value # y coordinate of the top 
		bottom_ypos = row[2].value # y coordinate of the bottom of the region

		# Read top(left) and bottom(right) coordinate for heatmap
		if i == 0:
			top = top_ypos # y coordinate of the top end of the embryo	
			bottom = bottom_ypos
		else: # bottom, replace previous assignment
			bottom = bottom_ypos # y coordinate of the bottom end of the embryo		
			
		# Plot slice boundaries
		for j in range(4,file_len_slice):
			row = list(worksheet_slice.row(j))
			top_left_xpos = row[3].value
			top_right_xpos = row[4].value
			bottom_left_xpos = row[1].value
			bottom_right_xpos = row[2].value
			
			if j == 4:	# get left most x-pos from first slice
				bottom_leftmost_xpos = bottom_left_xpos
				top_leftmost_xpos = top_left_xpos
			
			if j == file_len_slice-1:	# get right most x-pos from last slice
				bottom_rightmost_xpos = bottom_right_xpos
				top_rightmost_xpos = top_right_xpos

			# Vertical lines		
			ax_her1.plot([bottom_right_xpos, top_right_xpos], [bottom_ypos, top_ypos], c="k", alpha=1, linewidth=lw)
			ax_her1.plot([bottom_left_xpos, top_left_xpos], [bottom_ypos, top_ypos], c="k", alpha=1, linewidth=lw)		
			ax_her7.plot([bottom_right_xpos, top_right_xpos], [bottom_ypos, top_ypos], c="k", alpha=1, linewidth=lw)
			ax_her7.plot([bottom_left_xpos, top_left_xpos], [bottom_ypos, top_ypos], c="k", alpha=1, linewidth=lw)

		
		left = min(left,top_leftmost_xpos,bottom_leftmost_xpos) # x coordinate of the left end of the embryo
		right = max(right,top_rightmost_xpos,bottom_rightmost_xpos) # x coordinate of the right end of the embryo
		
		# Label left/right
		if i == 0: # left
			ax_her1.text((top_leftmost_xpos+top_rightmost_xpos)/2, top_ypos+20, "L")
			ax_her7.text((top_leftmost_xpos+top_rightmost_xpos)/2, top_ypos+20, "L")
			if len(worksheets) <= 1:
				ax_her1.text((bottom_leftmost_xpos+bottom_rightmost_xpos)/2, bottom_ypos-30, "R")
				ax_her7.text((bottom_leftmost_xpos+bottom_rightmost_xpos)/2, bottom_ypos-30, "R")
			
		else: # right, replace previous assignment for bottom value
			ax_her1.text((bottom_leftmost_xpos+bottom_rightmost_xpos)/2, bottom_ypos-30, "R")
			ax_her7.text((bottom_leftmost_xpos+bottom_rightmost_xpos)/2, bottom_ypos-30, "R")
			
		# Horizontal lines
		ax_her1.plot([bottom_leftmost_xpos, bottom_rightmost_xpos], [bottom_ypos, bottom_ypos], c="k", alpha=1, linewidth=lw)
		ax_her1.plot([top_leftmost_xpos, top_rightmost_xpos], [top_ypos, top_ypos], c="k", alpha=1, linewidth=lw)
		ax_her7.plot([bottom_leftmost_xpos, bottom_rightmost_xpos], [bottom_ypos, bottom_ypos], c="k", alpha=1, linewidth=lw)
		ax_her7.plot([top_leftmost_xpos, top_rightmost_xpos], [top_ypos, top_ypos], c="k", alpha=1, linewidth=lw)
		
		# Read individual cell data
		worksheet = worksheets[i]
		file_len =  worksheet.nrows
		
		cell_xpos = []
		cell_ypos = []
		her1_level = []
		her7_level = []
		for j in range(4,file_len):
			row = list(worksheet.row(j))
			cell_xpos.append(row[0].value)
			cell_ypos.append(row[1].value)
			her1_level.append(row[2].value) 
			her7_level.append(row[3].value)


		threshold_her1 = (max(her1_level) - min(her1_level))*0.3 + min(her1_level)
		threshold_her7 = (max(her7_level) - min(her7_level))*0.3 + min(her7_level)
		
		# Plot low expression cells first
		for j in range(len(cell_xpos)):
			if her1_level[j]<threshold_her1:			
				ax_her1.scatter(cell_xpos[j], cell_ypos[j], s=cs, c="black", alpha=0.2, edgecolors='none')		
			if her7_level[j]<threshold_her7:
				ax_her7.scatter(cell_xpos[j], cell_ypos[j], s=cs, c="black", alpha=0.2, edgecolors='none')		
		
		# Plot high expression cells on top of low expression cells
		for j in range(len(cell_xpos)):
			if her1_level[j]>=threshold_her1:			
				ax_her1.scatter(cell_xpos[j], cell_ypos[j], s=cs, c="b", alpha=1, linewidth=0.2)	
			if her7_level[j]>=threshold_her7:
				ax_her7.scatter(cell_xpos[j], cell_ypos[j], s=cs, c="r", alpha=1, linewidth=0.2)
		"""
		# mean threshold
		threshold_her7 = (max(her7_level) - min(her7_level))/10
		threshold_her1 = (max(her1_level) - min(her1_level))/10

		for j in range(len(cell_xpos)):
			for i in range(10):
				if her1_level[j]<= min(her1_level)+(i+1)*threshold_her1:
					ax_her1.scatter(cell_xpos[j], cell_ypos[j], s=cs, c=colorblue[i], alpha=0.8, edgecolors='none')
					break
			
			for i in range(10):
				if her7_level[j]<= min(her7_level)+(i+1)*threshold_her7:
					ax_her7.scatter(cell_xpos[j], cell_ypos[j], s=cs, c=colorred[i], alpha=0.8, edgecolors='none')
					break

		#by rank threshold
		sorted_her7 = sorted(her7_level)
		sorted_her1 = sorted(her1_level)
		step = int(len(sorted_her7)/10)
		threshold_her7 = []
		threshold_her1 = []
		for i in range(9):
			threshold_her7.append(sorted_her7[step*(i+1)])
			threshold_her1.append(sorted_her1[step*(i+1)])
		threshold_her7.append(sorted_her7[-1])
		threshold_her1.append(sorted_her1[-1])

		for j in range(len(cell_xpos)):
			for i in range(10):
				if her1_level[j]<= threshold_her1[i]:
					ax_her1.scatter(cell_xpos[j], cell_ypos[j], s=cs, c=colorblue[i], alpha=1, edgecolors='none')
					break

			for i in range(10):
				if her7_level[j]<= threshold_her7[i]:
					ax_her7.scatter(cell_xpos[j], cell_ypos[j], s=cs, c=colorred[i], alpha=1, edgecolors='none')
					break
		"""
		
	ax_her1.set_xlim(-50+left,right+50)
	ax_her7.set_xlim(-50+left,right+50)
	ax_her1.set_ylim(-50+bottom,top+50)
	ax_her7.set_ylim(-50+bottom,top+50)
	ax_her1.set_aspect('equal') # x and y axes have the same scaling
	ax_her7.set_aspect('equal') # x and y axes have the same scaling
	ax_her1.set_axis_off() # no axes
	ax_her7.set_axis_off() # no axes
	fig_her1.subplots_adjust(left=0, bottom=0, right=1, top=1,  wspace=0, hspace=0)	# no margin
	fig_her7.subplots_adjust(left=0, bottom=0, right=1, top=1,  wspace=0, hspace=0)	# no margin
	fig_her1.savefig(directory + "/heatmap_her1.png", format = "png", dpi=300)
	fig_her7.savefig(directory + "/heatmap_her7.png", format = "png", dpi=300)

def usage():
	print("create_heatmap.py: Invalid command-line arguments.")
	print("Format: python create_heatmap.py -i <input Excel file:cell_info> -s <input Excel file:slice info> -d <output directory>")
	print("Example: python create_heatmap.py -i ../wildtype/output/embryo1/cells.xls -s ../wildtype/output/embryo1/sliceInfo.xls -d ../wildtype/output/embryo1")
	exit(1)

main()
