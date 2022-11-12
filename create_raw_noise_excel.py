"""
Create an Excel file containing each slice's noise data for all embryos in a given genetic background
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
from xlrd import XLRDError

def main():
	# Check input
	if not shared.isInt(sys.argv[1]):
		print('create_raw_noise_excel.py: Number of embryos must be an integer.')
		exit(1)
	elif int(sys.argv[1])<=0:
		print('create_raw_noise_excel.py: Number of embryos must be larger than zero.')
		exit(1)		
	
	num_embryos = int(sys.argv[1])
	    	        
	if len(sys.argv)==num_embryos+3:
		inputs = sys.argv[2:num_embryos+2]
		directory = sys.argv[num_embryos+2]		
	else:
		usage()
		
	shared.ensureDir(directory)
	
	write_workbook = xlwt.Workbook(encoding="ascii")
	
	slice_mean_her1 = []
	slice_mean_her7 = []	
	slice_mean_her = []
	slice_mean_hm = []			
	slice_intrinsic_noise = []
	slice_extrinsic_noise = []
	slice_total_noise = []
	
	labels = ["Slice #", "Her1 mean", "Her7 mean", "Her mean", "Harmonic mean", "Intrinsic noise", "Extrinsic noise", "Total noise"]	
	for i in range(len(inputs)): # embryo
		# Open embryo data
		if not os.path.isfile(inputs[i]):
			print('create_raw_noise_excel.py: File "'+inputs[i]+'" does not exist.')
			exit(1)		
		try:
			read_workbook = xlrd.open_workbook(inputs[i],'r')
		except XLRDError as e:
			print('create_raw_noise_excel.py: Cannot open file "'+inputs[i]+'".')
			exit(1)	
		read_worksheets = read_workbook.sheets()
				
		write_worksheet = write_workbook.add_sheet(str(i+1))	
		write_worksheet.write(0,0,"Left")
		write_worksheet.write(0,9,"Right")
		for j in range(len(labels)):
			write_worksheet.write(1,j,labels[j])
			write_worksheet.write(1,j+9,labels[j])
					
		for j in range(len(read_worksheets)): # region 
			read_worksheet = read_worksheets[j]
			file_len =  read_worksheet.nrows	
						
			for k in range(1,file_len): # slice
				row = list(read_worksheet.row(k))			
				
				# Blank as default values in case data is invalid
				her1_mean = ""
				her7_mean = ""
				her_mean = ""
				harmonic_mean = ""
				intrinsic_noise = ""
				extrinsic_noise = ""
				total_noise = ""				
					
				if isinstance(row[1].value,float): # valid slice
					num_cells = int(row[1].value) # number of cells within this slice											
					her1 = []
					her7 = []
					
					for l in range(num_cells): # background subtraction
						# Take the cell's data only if its expression levels are positive after background subtraction	
						if row[8+2*l].value>0 and row[8+2*l+1].value>0:
							her1.append(row[8+2*l].value)
							her7.append(row[8+2*l+1].value)
										
					if len(her1)>=3: # valid only if more than 3 cells											
						her1_mean = numpy.mean(her1)
						her7_mean = numpy.mean(her7)
						her_mean = her1_mean + her7_mean
						harmonic_mean = 2/(1/her1_mean+1/her7_mean)
																		
						slice_mean_her1.append(her1_mean) # her1 mean
						slice_mean_her7.append(her7_mean) # her7 mean
						slice_mean_her.append(her_mean) # her mean (her1 + her7)
						slice_mean_hm.append(harmonic_mean) # harmonic mean of her1 and her7
						
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
						
						# Total noise
						total_noise = intrinsic_noise+extrinsic_noise
						slice_total_noise.append(total_noise)
				
				# Write slice data to Excel
				line = [k, her1_mean, her7_mean, her_mean, harmonic_mean, intrinsic_noise, extrinsic_noise, total_noise]
				for l in range(len(line)):					
					write_worksheet.write(k+1,l,line[l]) if j == 0 else write_worksheet.write(k+1,l+9,line[l]) 
	
	# Write combined data														
	write_worksheet = write_workbook.add_sheet("Combined")	
	for i in range(len(labels)):
		write_worksheet.write(0,i,labels[i])
	for i in range(len(slice_mean_her1)):
		line = [i+1, slice_mean_her1[i], slice_mean_her7[i], slice_mean_her[i], slice_mean_hm[i], slice_intrinsic_noise[i], slice_extrinsic_noise[i], slice_total_noise[i]]
		for j in range(len(line)):
			write_worksheet.write(i+1,j,line[j])		
	write_workbook.save(directory + "/raw_noise.xls")
	
def usage():
	print("create_raw_noise_excel.py: Invalid command-line arguments")
	print("Format: python create_raw_noise_excel.py <number of embryos> <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls> <output directory>")
	print("Example: python create_raw_noise_excel.py 20 ../wildtypefulldataset/output/embryo1/slices.xls ../wildtypefulldataset/output/embryo2/slices.xls ... ../wildtypefulldataset/output/embryo20/slices.xls ../wildtypefulldataset/output")
	exit(1)
	
main()
