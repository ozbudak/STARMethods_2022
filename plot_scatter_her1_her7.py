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

import sys, shared, os
import numpy, math
import matplotlib.pyplot as plt
import matplotlib.legend_handler as handler
import matplotlib.lines as mlines
import xlrd, xlwt
from xlrd import XLRDError
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [1,2,5,10,20,30,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def get_her_bounds(ws_bounds, num_bins):
	bounds = []
	for i in range(num_bins):
		bound = shared.toFlo(ws_bounds.cell_value(1, i))
		bounds.append(bound)
	return bounds
	
def get_row_data(row, start_index, end_index):
	results = []
	i = start_index
	while i < end_index: 
		results.append(shared.toFlo(row[i].value))
		i += 1
	return results
	
def plot_bin_scatter_her1_her7(ws_data, num_bins, outdir, bounds):
	# 1. Create num_bins figure objects
	figs = []
	axes = []
	for i in range(num_bins):
		figure = plt.figure(figsize = (6,6), dpi = 300)
		ax = figure.add_subplot(111)
		figs.append(figure)
		axes.append(ax)
	# 2. Plot for each slice from all embryos
	num_rows = ws_data.nrows
	for i in range(1, num_rows):
		row = list(ws_data.row(i))
		avg_her = shared.toFlo(row[0].value)
		for j in range(num_bins):
			if (avg_her <= bounds[j]): # this slice's avg_her is <= the upper bound of the bin --> found the correct bin
				num_cells = shared.toInt(row[7].value)
				her1 = get_row_data(row, 8, (8 + num_cells))
				her7 = get_row_data(row, (8 + num_cells), (8 + 2 * num_cells))
				axes[j].scatter(her1, her7, alpha = 0.1, color = '#000000')
				break
	# 3. Fix the ticks for axises of figures, and create axis labels
	tick_interval = determineTickInterval(200, 5)
	ticks = numpy.arange(0, 200, tick_interval)
	for i in range(num_bins):
		axes[i].set_xlim([0,200])
		axes[i].set_ylim([0,200])
		axes[i].set_yticks(ticks)
		axes[i].set_xticks(ticks)
		axes[i].tick_params(direction='in')
		axes[i].set_xlabel("$\mathit{her1}$ mRNA")
		axes[i].set_ylabel("$\mathit{her7}$ mRNA")
		(figs[i]).subplots_adjust(left=0.1, bottom=0.1, right=.95, top=.95,  wspace=None, hspace=None)

	# 4. Save figures for all bins
	for i in range (num_bins): 
		(figs[i]).savefig(outdir + "/her1_her7_bin" + str(i) + ".png", format = "png", dpi = 300)

def plot_all_scatter_her1_her7(ws_data, outdir):
	fig = plt.figure(figsize= (6,6), dpi = 300)
	ax = fig.add_subplot(111)
	num_rows = ws_data.nrows
	for i in range(1, num_rows):
		row = list(ws_data.row(i))
		num_cells = shared.toInt(row[7].value)
		her1 = get_row_data(row, 8, (8 + num_cells))
		her7 = get_row_data(row, (8 + num_cells), (8 + num_cells * 2))
		ax.scatter(her1, her7, alpha = 0.1, color = '#000000')
	tick_interval = determineTickInterval(200, 5)
	ticks = numpy.arange(0, 200, tick_interval)
	ax.set_xlim([0,200])
	ax.set_ylim([0,200])
	ax.set_xticks(ticks)
	ax.set_yticks(ticks)
	ax.tick_params(direction='in')
	ax.set_xlabel("$\mathit{her1}$ mRNA")
	ax.set_ylabel("$\mathit{her7}$ mRNA")
	fig.subplots_adjust(left=0.1, bottom=0.1, right=.95, top=.95,  wspace=None, hspace=None)
	fig.savefig(outdir + "/her1_her7_all.png", format = "png", dpi = 300)
	
def main():
	if len(sys.argv) != 4:
		usage()
	# Open the combined_slices excel file
	if not os.path.isfile(sys.argv[1]):
		print("plot_scatter_her1_her7.py: File " + sys.argv[1] + " does not exist.")
		exit(1)
	try:
		wb = xlrd.open_workbook(sys.argv[1], 'r')
	except XLRDError as e:
		print("plot_scatter_her1_her7.py: Cannot open file " + sys.argv[1] + " =.=!")
		exit(1)
	num_bins = shared.toInt(sys.argv[2])
	outdir = sys.argv[3]
	shared.ensureDir(outdir)
	
	# open worksheets
	ws_data = wb.sheet_by_name("aggregate_data")
	ws_bounds = wb.sheet_by_name("bin_bounds")
	# get the upper bounds for average her concentration of all bins
	bounds = get_her_bounds(ws_bounds, num_bins)
	# plot the slices into different figures, corresponding to different bins
	plot_bin_scatter_her1_her7(ws_data, num_bins, outdir, bounds)
	# scatter plot of her1 and her7 in all bins
	plot_all_scatter_her1_her7(ws_data, outdir)
	
def usage():
	print("Plot_scatter_her1_her7.py: Plotting the concentrations of her1 and her7 mRNAs in cells across all embryos")
	print("python plot_scatter_her1_her7.py <inputFile: combined_slices.xls-- result of plot_noise.py> <number of bins> <output directory to store the plots>")
	print("Example: python plot_scatter_her1_her7.py ../wildtypefulldataset/output/combined_slices.xls 5 ../wildtypefulldataset/output")
	exit(1)

main()
