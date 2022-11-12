"""
Create Excel files containing each slice's RNA expression levels 
after background noise subtraction for all embryos in a given genetic background
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
import numpy
import xlrd, xlwt
import matplotlib.pyplot as plt


def plot_her(combine_her1, combine_her7, combine_total_her, directory):
	plt.subplot(211)
	#plt.hist(combine_her7, 50, normed=1, facecolor='lightblue', alpha=1.0, label='her7')
	#plt.hist(combine_her1, 50, normed=1, facecolor='orange', alpha=0.5, label='her1')
	plt.hist(combine_her7, 50, facecolor='lightblue', alpha=1.0, label='her7')
	plt.hist(combine_her1, 50, facecolor='orange', alpha=0.5, label='her1')
	plt.legend()
	plt.tick_params(direction='in')
	plt.ylabel('Frequency')
	#plt.xlabel('mRNA expression level')
	#plt.savefig(directory + "/combineher1her7hist.png", format = "png")

	# Plot total her expression distribution
	plt.subplot(212)
	#plt.hist(combine_total_her, 50, normed=1, facecolor='orange', alpha=1.0, label='Total her')	
	plt.hist(combine_total_her, 50, facecolor='orange', alpha=1.0, label='Total her')	
	plt.legend()
	plt.tick_params(direction='in')
	plt.ylabel('Frequency')
	plt.xlabel('mRNA expression level')
	plt.savefig(directory + "/combinetotalherhist.png", format = "png", dpi=300)
	
def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print('create_raw_expression_excel.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print('create_raw_expression_excel.py: Number of embryos must be larger than zero.')
		exit(1)		
	
	num_embryos = int(sys.argv[1])
	    	        
	if len(sys.argv)==num_embryos+3:
		inputs = sys.argv[2:num_embryos+2]
		directory = sys.argv[num_embryos+2]		
	else:
		usage()
		        
	shared.ensureDir(directory) 	
	
	# Set up workbooks to write
	write_workbook_after = xlwt.Workbook(encoding="ascii")

	# Collect all her1, her7 mRNA count for histogram plotting
	combine_her1 = []
	combine_her7 = []
	combine_total_her = []
		
	for i in range(len(inputs)): # embryo
		# Open embryo data
		if not os.path.isfile(inputs[i]):
			print('create_raw_expression_excel.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			read_workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print('create_raw_expression_excel.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		read_worksheets = read_workbook.sheets()
		
		# Set up worksheets to write		
		write_worksheet_after = write_workbook_after.add_sheet(str(i+1))
		labels = ["Region","Slice #","Gene"]
		for j in range(len(labels)):
			write_worksheet_after.write(0,j,labels[j])		
		row_num = 1	
		
		for j in range(len(read_worksheets)): # region 
			read_worksheet = read_worksheets[j]
			file_len =  read_worksheet.nrows	
			region = "L" if j == 0 else "R"	
				
			for k in range(1,file_len): # slice
				row = list(read_worksheet.row(k))
				her1_after = []
				her7_after = []
				total_her_after = []
					
				if isinstance(row[1].value,float): # valid slice
					num_cells = int(row[1].value) # number of cells within this slice
					
					for l in range(num_cells):				
						her1_after.append(row[8+2*l].value) 
						her7_after.append(row[8+2*l+1].value)
						total_her_after.append(row[8+2*l].value + row[8+2*l+1].value)
				
				if len(her1_after)==0:
					her1_after = ["Too few cells"]
					her7_after = ["Too few cells"]	
				
				line_her1_after = [region, k, "Her1"] + her1_after
				line_her7_after = [region, k, "Her7"] + her7_after
				
				for l in range(len(line_her1_after)):
					write_worksheet_after.write(row_num,l,line_her1_after[l])	
					write_worksheet_after.write(row_num+1,l,line_her7_after[l])
				row_num+=2 # used two rows for her1 and her7	

				combine_her1 = combine_her1 + her1_after
				combine_her7 = combine_her7 + her7_after
				combine_total_her = combine_total_her + total_her_after
					
	write_workbook_after.save(directory + "/raw_expression_afterbackgroundsub.xls")
	combine_her1 = [i for i in combine_her1 if isinstance(i, (int, float))]
	combine_her7 = [i for i in combine_her7 if isinstance(i, (int, float))]
	plot_her(combine_her1, combine_her7, combine_total_her, directory)
	
def usage():
	print("create_raw_expression_excel.py: Invalid command-line arguments")
	print("Format: python create_raw_expression_excel.py <number of embryos> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print("Example: python create_raw_expression_excel.py 20 ../wildtypefulldataset/output/embryo1/slices.xls ../wildtypefulldataset/output/embryo2/slices.xls ... ../wildtypefulldataset/output/embryo20/slices.xls ../wildtypefulldataset/output")
	exit(1)

main()
