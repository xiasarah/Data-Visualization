#!/usr/bin/env python
# coding: utf-8

# In[87]:


import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 
import matplotlib.patheffects as PathEffects

plt.style.use('BME163.mplstyle')

parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')
parser.add_argument('--inFile1','-c',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile2','-p',type = str,action = 'store',help = 'input file')

args = parser.parse_args()
outFile = args.outFile
inFile1 = args.inFile1
inFile2 = args.inFile2


figureWidth = 8
figureHeight = 4

panel1Width = 2
panel1Height = 2

plt.figure(figsize=(figureWidth,figureHeight))

panel1 = plt.axes([0.1,0.2,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical

panel1.set_ylim(-40,30)
panel1.set_xlim(-30,30)
panel1.set_xticks([-20,0,20])

plt.xlabel("tSNE 2")
plt.ylabel("tSNE 1")




# reading files


cellDict={}

for line in open(inFile1):
    a = line.rstrip().split('\t')
    cell = a[0]
    cellType = a[1]
    barcode = a[2]
    cellDict[barcode] = cellType

    

xList_monocyte = []
yList_monocyte = []

xList_tCell = []
yList_tCell = []

xList_bCell = []
yList_bCell = []


for line in open(inFile2):
    a = line.rstrip().split()
    seq = a[0]
    x = float(a[1])
    y = float(a[2])
    
    if cellDict.get(seq) == 'monocyte':
        xList_monocyte.append(x)
        yList_monocyte.append(y)
        
    elif cellDict.get(seq) == 'tCell':
        xList_tCell.append(x)
        yList_tCell.append(y)
    
    elif cellDict.get(seq) == 'bCell':
        xList_bCell.append(x)
        yList_bCell.append(y)

# plotting

panel1.plot(xList_bCell,yList_bCell,
                marker = 'o',
                markersize = 4,
                linewidth = 0,
                linestyle = '--',
                markeredgewidth = 0,
                markeredgecolor = (70/250,130/250,180/250),
                markerfacecolor = (70/250,130/250,180/250)
                )

panel1.plot(xList_tCell,yList_tCell,
                marker = 'o',
                markersize = 4,
                linewidth = 0,
                linestyle = '--',
                markeredgewidth = 0,
                markeredgecolor = (221/250,160/250,221/250),
                markerfacecolor = (221/250,160/250,221/250)
                )

panel1.plot(xList_monocyte,yList_monocyte,
                marker = 'o',
                markersize = 4,
                linewidth = 0,
                linestyle = '--',
                markeredgewidth = 0,
                markeredgecolor = (102/250,205/250,170/250),
                markerfacecolor = (102/250,205/250,170/250)
                )

# panel text

x_monocyte = np.median(xList_monocyte)
y_monocyte = np.median(yList_monocyte)

x_tCell = np.median(xList_tCell)
y_tCell = np.median(yList_tCell)

x_bCell = np.median(xList_bCell)
y_bCell = np.median(yList_bCell)


monocyteText = panel1.text(x_monocyte,y_monocyte,'monocyte',fontsize=8,va='center',ha='center')
tCellText = panel1.text(x_tCell,y_tCell,'tCell',fontsize=8,va='center',ha='center')
bCellText = panel1.text(x_bCell,y_bCell,'bCell',fontsize=8,va='center',ha='center')

monocyteText.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
tCellText.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])
bCellText.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])





plt.savefig(outFile, dpi = 600)


# In[ ]:




