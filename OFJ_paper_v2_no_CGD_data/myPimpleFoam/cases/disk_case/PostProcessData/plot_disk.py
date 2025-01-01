# -*- coding: utf-8 -*-
"""
@author: Authors of manuscript #131 submitted to the OpenFOAM Journal.
 Author names were suppressed to allow the double-blind review process, they will be added at the end of the review
"""

#%% 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.pyplot import cm
import sys as sys

#%%
fileA2U    = sys.argv[1]
fileSol    = sys.argv[2]
listOfCmpt = [str(i) for i in sys.argv[3].split(',')]
listOfCmpt.sort()
saveFigPath='/'.join(fileA2U.split('/')[0:-1])

# Dictionary to map key, values
changeDict = {'0':'A_{11}', 
              '1':'A_{12}',
              '2':'A_{13}',
              '3':'A_{22}',
              '4':'A_{23}',
              '5':'A_{33}'}

changeDictOF = {'A_{11}': 'A2_0', 
              'A_{12}': 'A2_1',
              'A_{13}': 'A2_2',
              'A_{22}': 'A2_3',
              'A_{23}': 'A2_4',
              'A_{33}': 'A2_5'}

legendDict = {'A_{11}': 'A_{xx}', 
              'A_{12}': 'A_{xy}',
              'A_{13}': 'A_{xz}',
              'A_{22}': 'A_{yy}',
              'A_{23}': 'A_{yz}',
              'A_{33}': 'A_{zz}'}

# Read data from user
dataOF = pd.read_csv(fileA2U)  
dataPython = pd.read_csv(fileSol)

# round z values in OpenFOAM some cells have minor values and pandas will not split this correct. 9 decimal cases should be enough
dataOF['z'] = dataOF['z'].round(9)

# Center data from OF to half thickness
h = 0.0015 # half thickness
dataOF['z_normalized'] = (dataOF['z'] - h)/h

# Filter data with positive values of thickness
dataOF['filter'] = dataOF['z_normalized'].apply(lambda x: True if x>=0 else False)
dataOF = dataOF.loc[dataOF['filter'], :]

# Apply the same filter to Python data
dataPython['z_normalized'] = dataPython['z']/h # Data from python is already centered
dataPython['filter'] = dataPython['z_normalized'].apply(lambda x: True if x>=0 else False)
dataPython = dataPython.loc[dataPython['filter'], :]


# Split data by unique z values. Sort it in ascending manner by the radius
outOF = [y.sort_values('x') for x, y in dataOF.groupby('z_normalized')]
outPython = [y.sort_values('r') for x, y in dataPython.groupby('z_normalized')]

# Check if the data same is correct
shapeToCheck = outOF[0].shape
for i in outOF:
    if (i.shape != shapeToCheck):
        print("something is wrong with the data")
        sys.exit(0)

if (len(outPython) != len(outOF)):
    raise ValueError('Something is wrong with the data')

# Number of lines to appear in plot- Find the closest values from 0.1-0.9 
targetValues = np.arange(0.1, 1, 0.1)
indexList = np.zeros(len(targetValues))

for idx, value in enumerate(targetValues):
    indexList[idx] = min(range(len(outOF)), key=lambda i: max(abs(outOF[i]['z_normalized'] - value)))

indexList=indexList.astype(int)

# Number of points to display in plot
numberOfPointsInPlot = 128
indexListPlot = np.round(np.linspace(0, shapeToCheck[0] - 1, numberOfPointsInPlot )).astype(int)


# Plot settings
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams['axes.linewidth'] = 2
fZ=14

colors = cm.rainbow(np.linspace(0, 1, len(targetValues)))
colorID=0

for cmpt in range(len(listOfCmpt)):
    
    plt.figure()
    
    A2cmpt = changeDict[listOfCmpt[cmpt]]
    A2legend = legendDict[A2cmpt]
    
    # Will make two loops for the legend to the correctly formated
    for idx in indexList:
        c = colors[colorID]
        plt.plot(
                 outPython[idx]['r'], 
                 outPython[idx][A2cmpt],
                 color=c,
                 path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()],
                 label='python',
                 zorder=0)
        colorID += 1
    colorID=0
    
    for idx in indexList:
        c = colors[colorID]
        plt.plot(
                    outOF[idx]['x'].iloc[indexListPlot],
                    outOF[idx][changeDictOF[A2cmpt]].iloc[indexListPlot],
                    linestyle='None',
                    color=c,
                    marker='o',
                    markersize=4,
                    markeredgecolor='black',
                    label='z/b= ' + str(round(outOF[idx]['z_normalized'].iloc[0], 3)),
                    zorder=1
                    )
        colorID += 1
            
        # Dummy check
    if(np.max(np.abs(outPython[idx]['z_normalized'] - outOF[idx]['z_normalized']))>1e-10):
            raise ValueError('thicknesses are different')
    colorID=0
    
    # Additional plot settings
    plt.xlim([0.01, 0.12])
    plt.legend(loc='lower right', ncol=2, fontsize=fZ-2)
    plt.xlabel("r [m]", fontsize=fZ)
    plt.ylabel(rf'${A2legend}$', fontsize=fZ)
    plt.yticks(fontsize = fZ+2)
    plt.xticks(fontsize = fZ+2)
    plt.grid(visible=True, which='major', color='black', linestyle='-')
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='darkgrey', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(saveFigPath +'/' + A2cmpt +".svg")


