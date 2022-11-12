"""
Run python scripts to create figures that compare data from different genetic backgrounds
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

def checkFile(f): # check if a given file exists and exit program if not
	if not os.path.isfile(f):
		print ("compare_geneticbackgrounds.py: File '"+f+"' does not exist.")
		exit(1)

def main():
	# Check input
	try:
		if not shared.isInt(sys.argv[1]):
			print ('compare_geneticbackgrounds.py: Number of genetic backgrounds must be an integer.')
			exit(1)
		elif not shared.isInt(sys.argv[2]):
			print ('compare_geneticbackgrounds.py: Number of bins for noise plots must be an integer.')
			exit(1)

		num_geneticbackgrounds = int(sys.argv[1]) # number of genetic backgrounds to combine (currently 3, can be increased)
		num_bins = int(sys.argv[2])
		output_directory = sys.argv[3]
	except: 
		usage()
		
	if len(sys.argv)==3*num_geneticbackgrounds+4:
		folders = sys.argv[4:num_geneticbackgrounds+4] # folders where the input files are located 
		for folder in folders:
			if not os.path.isdir(folder):
				print ("compare_geneticbackgrounds.py: Folder '"+folder+"' does not exist.")
				exit(1)
		names = sys.argv[num_geneticbackgrounds+4:2*num_geneticbackgrounds+4] # strings containing genetic background names
		colors = sys.argv[2*num_geneticbackgrounds+4:] # colors for plotting genetic backgrounds
	else:
		usage()	
	
	
	### Compare her1 and her7 spatial amplitudes ###
	print ("Plotting spatial amplitudes from different genetic backgrounds...")
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/spatial_amplitude.xls"		
		checkFile(f)
		files.append(f)
	command = ["python","compare_spatial_amplitude.py",str(num_geneticbackgrounds), output_directory] + files + names + colors
	if 1==call(command):
		exit(1)
	
	
	
    ###Compare burst size and frequencies ###
	'''print ("Plotting estimated burst parameters from different genetic backgrounds...")
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/estimated_burst_parameters.xls"
		checkFile(f)
		files.append(f)
	command = ["python","compare_burst_parameters.py",str(num_geneticbackgrounds), output_directory] + files + names + colors
	if 1==call(command):
		exit(1)'''
	
	### Compare the Fano factor ###
	print ("Plotting the Fano factor from different genetic backgrounds...")
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/fano_factor.xls"
		checkFile(f)
		files.append(f)
	command = ["python","compare_fano_factor.py", str(num_geneticbackgrounds), output_directory] + files + names + colors
	if 1==call(command):
		exit(1)
	
	### Compare noise levels ###
	print ("Plotting noise from different genetic backgrounds...")
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/noise.xls"
		checkFile(f)
		files.append(f)
	command = ["python","compare_noise.py", str(num_geneticbackgrounds), '5', output_directory] + files + names + colors 
	if 1==call(command):
		exit(1)		
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/raw_noise.xls"
		checkFile(f)
		files.append(f)
	# Noise is normalized to the first dataset
	command = ["python","compare_noise_bar.py", str(num_geneticbackgrounds), output_directory] + files + names + colors + [files[0]]
	if 1==call(command): 
		exit(1)
	
	### Compare CV^2 ###
	print ("Plotting CV^2 from different genetic backgrounds...")
	files = []
	for i in range(num_geneticbackgrounds):
		f = folders[i]+"/CVsquared.xls"
		checkFile(f)
		files.append(f)
	command = ["python","compare_grouped_CVsquared.py", str(num_geneticbackgrounds), output_directory] + files + names + colors
	if 1==call(command):
		exit(1)	
	
	command = ["python","compare_CVsquared.py", str(num_geneticbackgrounds), output_directory] + files + names + colors
	if 1==call(command):
		exit(1)	

	command = ["python","plot_CVsquare_fillArea.py", output_directory] + colors
	if 1==call(command):
		exit(1)

	command = ["python","compare_CVsquared_fillArea.py", output_directory] + colors
	if 1==call(command):
		exit(1)

def usage():
	print ("compare_geneticbackgrounds.py: Invalid command-line arguments.")
	print ("Format: python COMPARE_GENETICBACKGROUNDS.py <number of genetic backgrounds> <number of bins for noise plots> <comparison outcome folder> <folders containing files to be compared (i.e. directory of “Output” folders of each genotype)> <genetic background names for labeling> <genetic background colors for plotting>")
	print ("Example: COMPARE_GENETICBACKGROUNDS.py 2 5 ../Compare/ ../wVol/b567xher17/output/ ../wVol/her17/output/ b567xher17 her17 b r")
	exit(1)
	
main()
