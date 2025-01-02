# -*- coding: utf-8 -*-
"""
@author: Authors of manuscript #131 published in the OpenFOAM Journal.
  doi: https://doi.org/10.51560/ofj.v5.131 

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import StrMethodFormatter

import sys as sys
import os

fileA2U    = sys.argv[1]
fileSol    = sys.argv[2]
listOfCmpt = [str(i) for i in sys.argv[3].split(',')]
listOfCmpt.sort()
plotTitle  = sys.argv[4]

saveFigPath='/'.join(fileA2U.split('/')[0:-1])



print(f"Using OF file: {os.path.abspath(fileA2U)}\n \
        Using python file: {fileSol }\n \
        With components: {listOfCmpt}\n ")

changeDict = {'0':'A_{xx}', 
              '1':'A_{xy}',
              '2':'A_{xz}',
              '3':'A_{yy}',
              '4':'A_{yz}',
              '5':'A_{zz}'}

dataA2U = pd.read_csv(fileA2U)  
dataSol = pd.read_csv(fileSol)

dataA2U['t_line'] = dataA2U["x"]/dataA2U["U_0"]

# Check where y value changes
dataA2U['diffY'] = (dataA2U['y'].diff() >= 1e-9)

# Split data by diffY
out = [dataA2U for _, dataA2U in dataA2U.groupby(dataA2U['diffY'].eq(1).cumsum())]

shapeToCheck = out[0].shape
for i in out:
    if (i.shape != shapeToCheck):
        print("something is wrong with the data")
        sys.exit(0)

# Compute cell centers for 32x32 mesh
nCells = 32
points= np.linspace(0, 1 , nCells+1)
uniformPoints=(points[1 : nCells+1] + points[0:nCells])*0.5
    
# Plot settings
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams['axes.linewidth'] = 2

fZ=17

limitsDict = {'A_{xx}':[0.3, 0.825], 
              'A_{xy}':[0, 0.25],
              'A_{xz}':[0, 1e-12],
              'A_{yy}':[0, 0.35],
              'A_{yz}':[0, 1e-12],
              'A_{zz}':[0, 0.35]}


for cmpt in range(len(listOfCmpt)):

    plt.figure()
    
    A2cmpt = 'A2'+'_'+ str(listOfCmpt[cmpt])
    A2Legend = changeDict[listOfCmpt[cmpt]]
    
    plt.plot(
                dataSol['time'], 
                dataSol[A2cmpt], 
                path_effects=[pe.Stroke(linewidth=5, foreground='black'), pe.Normal()],
                label=fr'python',
                linewidth=3,
                zorder=1
            )
    
    # Plots first 3 layers in the y-direction
    for idx, i in enumerate(out[:3]):
            indexListPlot = np.zeros(nCells)
            for cell in range(nCells):
                indexListPlot[cell]=np.abs(i['x'] - uniformPoints[cell]).argmin()
            
            plt.plot(
                        i['t_line'].iloc[indexListPlot],
                        i[A2cmpt].iloc[indexListPlot],
                        linestyle='None',
                        marker='o',
                        markersize=7,
                        markeredgecolor='black',
                        label=fr'y= '  + str(round(i['y'].iloc[0],4)) + ' [m]',
                        zorder=2
                        )

    plt.xlim([0, 63])  
    plt.ylim(limitsDict[A2Legend])
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    plt.legend(loc='lower right', fontsize=fZ-4)
    plt.xlabel("t' [s]", fontsize=fZ)
    plt.ylabel(rf'${ A2Legend }$', fontsize=fZ)
    plt.yticks(fontsize = fZ+2)
    plt.xticks(fontsize = fZ+2)
    plt.grid(visible=True, which='major', color='black', linestyle='-', zorder=0)
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='darkgrey', linestyle='--', alpha=0.3, zorder=0)
    plt.title(plotTitle, loc='right')
    plt.tight_layout()
    plt.savefig(saveFigPath +'/' + plotTitle + '_' + A2Legend + ".svg", transparent=True)

