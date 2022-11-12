"""
Create all figures and Excel files that combine data from all embryos in a given genetic background
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
from subprocess import call

DEFAULT_NUM_BIN = 5

def main():
	args = sys.argv[1:]
	num_args = len(args)
	req_args = [False]*3
	
	num_bins = DEFAULT_NUM_BIN
	if num_args >= 6:
		i = 0
		while i<num_args-1:
			option = args[i]
			value = args[i+1]			
			if (option == '-ne' or option == '--number-of-embryos') and shared.isInt(value):
				num_embryos = int(value)
				req_args[0] = True
				i+=2
			elif (option == '-nb' or option == '--number-of-bins') and shared.isInt(value):
				num_bins = int(value)
				i+=2
			elif option == '-d' or option == '--output-directory':
				directory = value
				req_args[1] = True
				i+=2
			elif req_args[0] and (option == '-i' or option == '--input-files') and ((num_args-7)==num_embryos):
				slice_files = args[i+1:i+1+num_embryos]
				for f in slice_files:
					if not os.path.isfile(f):
						print("combine_embryos.py: File "+f+" does not exist.")
						exit(1)
				req_args[2] = True	
				i+=num_embryos
			else:
				usage()
		for arg in req_args:
			if not arg:
				usage()
	else:
		usage()

	shared.ensureDir(directory)

	### Spatial amplitude ###
	print("Plotting spatial amplitude...")
	command = ["python","plot_spatial_amplitude.py",str(num_embryos)] + slice_files + [directory]
	if 1==call(command):
		exit(1)
	# (compare_spatial_amplitude.py can run after plot_spatial_amplitude.py is run for all genetic backgrounds)		

	'''### Burst size and frequency ###
	# 1. create_burst_data.py
	print("Creating data for estimate_burst_parameters.m...")
	command = ["python","create_burst_data.py",str(num_embryos)] + slice_files + [directory]
	if 1==call(command):
		exit(1)
	# 2. estimate_burst_parameters.m 
	print("Running estimate_burst_parameters.m on MATLAB...")
	command = ['/Applications/MATLAB_R2016a.app/bin/matlab','-nodesktop','-nosplash','-nodisplay','-r','estimate_burst_parameters(\''+directory+'/burst_data.xls\',\''+directory+'\')']
	if 1==call(command): # this will automatically open and run MATLAB
		exit(1)
	
	# 3. plot_estimated_burst_parameters.py using the output from estimate_burst_parameters.m 	
	print("Plotting estimated burst size and frequencies...")
	command = ["python","plot_estimated_burst_parameters.py",directory+"/burst_result.xls",directory]
	if 1==call(command):
		exit(1)'''	
	# (compare_burst_parameters.py can run after plot_estimated_burst_parameters.py is run for all genetic backgrounds)
	
	# Fano factor (to demonstrate burstiness) 
	command = ["python","plot_fano_factor.py",str(num_embryos)] + slice_files + [directory]
	print("Plotting fano factor...")
	if 1==call(command):
		exit(1)
	# (compare_fano_factor.py can run after plot_fano_factor.py is run for all genetic backgrounds)

	### Noise ###
	# Intrinsic and extrinsic noise
	print("Plotting intrinsic and extrinsic noise...")
	command = ["python","plot_noise.py",str(num_embryos), str(num_bins)] + slice_files + [directory]
	if 1==call(command):
		exit(1)	
	# (compare_noise.py can run after plot_noise.py is run for all genetic backgrounds)	
	
	### Scatter plot of her1 and her7 for all bins ####
	print("Plotting scatter plots of her1 vs her7 mRNAs in all bins ...")	
	command = ["python", "plot_scatter_her1_her7.py", directory + "/combined_slices.xls", str(num_bins), directory]
	if 1 == call(command):
		exit(1)
	
	# Spatial noise (coefficient of variation squared across space)	
	print("Plotting spatial noise (coefficient of variation squared across space)...")
	command = ["python","plot_CVsquared.py",str(num_embryos)] + slice_files + [directory]
	if 1==call(command):
		exit(1)	
	
	# (compare_grouped_CVsquared.py and compare_CV_squared.py can run after plot_CVsquared.py is run for all genetic backgrounds)		
	### Raw data Excel files ###
	command = ["python","create_raw_expression_excel.py",str(num_embryos)] + slice_files + [directory]
	print("Creating Excel files for RNA expression levels...")
	if 1==call(command):
		exit(1)	

	command = ["python","create_raw_spacial_noise_excel.py",str(num_embryos)] + slice_files + [directory]
	print("Creating Excel files for spacial noise...")
	if 1==call(command):
		exit(1)
	
	command = ["python","create_raw_noise_excel.py",str(num_embryos)] + slice_files + [directory]
	print("Creating Excel files for noise...")
	if 1==call(command):
		exit(1)

def usage():
	print("combine_embryos.py: Invalid command-line arguments.")
	print("Format: combine_embryos.py -ne <number of embryos> -nb <number of bins>  -d <output directory> -i <first embryo's slice.xls> <second embryo's slice.xls> ... <last embryo's slice.xls>")
	print("Example: python combine_embryos.py -ne 20 -d ../wildtypefulldataset/output -nb 5 -i ../wildtypefulldataset/output/embryo1/slices.xls \
	../wildtypefulldataset/output/embryo2/slices.xls .... ../wildtypefulldataset/output/embryo20/slices.xls")
	exit(1)
	
main()
