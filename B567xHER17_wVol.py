'''
Run python scripts to create all figures needed for wild-type embryos
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
import xlrd
from subprocess import call

############ THE FOLLOWING VALUES CAN BE CHANGED IF THE INPUT VALUES ARE CHANGED

# where input and output are located
folderIn = '../wVol/b567xher17/input'
folderOut = '../wVol/b567xher17/output'

## Reading embryo information from SampleInfo.xlsx and create  arrays - modified by leeyy
sampleInfo = folderIn + "/SampleInfo.xlsx" # named information file as "SampleInfo.xlsx" and put it at input folder

workbook = xlrd.open_workbook(sampleInfo,'r')
worksheet = workbook.sheet_by_name("Sheet1")
file_len = worksheet.nrows

num_embryos = file_len - 1
left_angles = [None] * num_embryos
right_angles = [None] * num_embryos
angle = 41.381 #angle of expression
delta_angle = 0.0403 #angle variation slope
CB = [None] * num_embryos
VARCB = [None] * num_embryos
YB = [None] * num_embryos
VARYB = [None] * num_embryos

for j in range(1, file_len): # skipped first line
    row = list(worksheet.row(j))    
    left_angles[j-1] = 180 + row[2].value   #column3: L angle
    right_angles[j-1] = row[3].value - 180  #column4: R angle
    CB[j-1] = row[4].value                  #column5: mean somite farred (HER1)
    VARCB[j-1] = row[6].value               #column7: somite variance farred (HER1)
    YB[j-1] = row[5].value                  #column6: mean somite red (HER7)
    VARYB[j-1] = row[7].value               #column8: somite variance red (HER7)
    
    
def main():
	commands = [] # list of commands for running embryo_analysis.py for each embryo
	slices_files = [] # list of slices.xls files from all embryos needed for combine_embryos.py
	# Putting commands into array commands, to be called later. If you change the name of input files, please modify
	# using the "-i" flag
	# Comment starting here if you want to skip embryo_analysis.py
	
	for i in range(1,num_embryos+1):
		commands.append(['python','embryo_analysis.py','-i',folderIn+'/b567xher17_wv_'+str(i)+'.xlsx','-d',folderOut+'/embryo'+str(i),'-a',str(angle),'-dA',str(delta_angle),'-n','2','-f','0','-m1',str(CB[i-1]),'-m7',str(YB[i-1])])
	# Process raw input data
	print('Analyzing b567xher17_wv embryos...')
	for command in commands:
		if 1==call(command):
			exit(1)
	
	# Comment ending here if you want to skip embryo_analysis.py
	# Plot heatmap and spatial expression for each embryo	
	print('Plotting heatmap and spatial expression for each embryo...')
	for i in range(1,num_embryos+1):
		# Comment starting here if you want to skip create_heatmap.py and plot_spatial_expression.py
		
		if 1==call(['python','create_heatmap.py','-i',folderOut+'/embryo'+str(i)+'/cells.xls','-d',folderOut+'/embryo'+str(i),'-s',folderOut+'/embryo'+str(i)+'/SliceInfo.xls']):
			exit(1)
		if 1==call(['python','plot_spatial_expression.py','-i',folderOut+'/embryo'+str(i)+'/slices.xls','-d',folderOut+'/embryo'+str(i)]):
			exit(1)
		
		# Comment ending here if you want to skip create_heatmap.py and plot_spatial_expression.py
	for i in range(1, num_embryos + 1):
		slices_files.append(folderOut+'/embryo'+str(i)+'/slices.xls')
	# Create figures that combine data from all embryos
	command = ['python','combine_embryos.py','-ne',str(num_embryos), '-nb', str(5),'-d',folderOut,'-i'] + slices_files
	if 1==call(command):
		exit(1)

	print('b567xher17_wv analysis... Done.')
main()
