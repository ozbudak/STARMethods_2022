# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 19:19:19 2018

@author: My
"""

import os, sys
import pandas as pd
import numpy, math
import matplotlib.pyplot as plt
from scipy.stats import sem
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [0.05,0.1,0.2,0.5,1,2]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 0.1

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=5)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)
    
def main():
    # Check input
    if len(sys.argv)>=3:
        output_directory = sys.argv[1]
        inputFile = output_directory + "/compare_CVsquared.xls"
        if not os.path.isfile(inputFile):
            print ("plot_CVsquared_fillArea.py: File '" + inputFile + "' does not exist.")
            exit(1)
        color = sys.argv[2:] # colors for plotting different genetic backgrounds
    else:
        usage()
                
    xls = pd.ExcelFile(inputFile)
    input_worksheet_names = ["Her1_plot_data", "Her7_plot_data", "Her_plot_data"]
    
    # Set up figure
    f = plt.figure(figsize=(6,12),dpi=600) # figure using only her on the x-axis 
    curPlot = 311
    
    for ws_name in input_worksheet_names:

        ax_her = f.add_subplot(curPlot)
        names = []
        worksheet = xls.parse(sheet_name=ws_name)
        geneNum = int((len(worksheet.columns)-1)/2)
        
        xmax = len(worksheet) # keep track of how long the x-axis should be
        #ymax = float('-inf') # keep track of how long the y-axis of each subplot should be
        #ymax = 1.1      # default maximum Y is determine manually

        legend_shapes = []
        legend_labels = []
        
        for i in range(geneNum):
            colName = worksheet.columns.values[i+1]
            names.append(colName[:-4])
            min_std = worksheet[colName]-worksheet[worksheet.columns.values[i+geneNum+1]]
            max_std = worksheet[colName]+worksheet[worksheet.columns.values[i+geneNum+1]]
            ax_her.plot(worksheet['Slide index'], worksheet[colName], color=color[i])
            ax_her.fill_between(worksheet['Slide index'], max_std, min_std, where=max_std>=min_std, facecolor=color[i], alpha=0.25)
            
            #ymax = max(ymax,numpy.array(max_std))
            legend_shapes.append(plt.Line2D(range(1), range(1), color="w", marker='o', markerfacecolor=color[i], markeredgecolor=color[i]))
            legend_labels.append(r''+names[i])
            
        #ymax = max(worksheet.columns.values)

        ax_her.legend(legend_shapes, legend_labels, numpoints=1, fontsize=12, loc=9, bbox_to_anchor=(0.5,1.1), ncol=geneNum, borderaxespad=0)
        ax_her.set_ylabel(r"$\textsf{CV}^{\textsf{\small{2}}}$ (%s)"%ws_name[:-10])
        ax_her.set_xlabel(r"Cell position (posterior - anterior)")	
        ax_her.xaxis.set_ticks(numpy.arange(0,xmax+1,10))
        ax_her.set_xlim(-1,xmax)
        #ax_her.yaxis.set_ticks(numpy.arange(0,ymax+1,determineTickInterval(ymax,5)))
        #ax_her.set_ylim(0,ymax*1.1)
        #updateTicklabels(ax_her)
        
        curPlot = curPlot + 1
        
    f.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=.95,  wspace=None, hspace=.3)	
    f.savefig(output_directory + "/CVsquare_fillArea.png", format = "png", dpi=600)
        
    
def usage():
	print ("plot_CVsquared_fillArea.py: Invalid command-line arguments")
	print ("Format: python plot_CVsquared_filled.py <output_directory> <list of colors to be used>")
	print ("Example: python compare_CVsquared.py 3 ../compare_output/WTdeltaCdeltaD ../wildtypefulldataset/output/CVsquared.xls ../deltacfulldataset/output/CVsquared.xls ../deltadfulldataset/output/CVsquared.xls Wildtype DeltaC DeltaD \#722AFF g r")
	exit(1)
    
main()
