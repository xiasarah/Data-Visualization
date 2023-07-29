#!/usr/bin/env python
# coding: utf-8

# In[75]:


import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 
import math

plt.style.use('BME163.mplstyle')

parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')
parser.add_argument('--inFile','-i',type = str,action = 'store',help = 'input file')

args = parser.parse_args()
outFile = args.outFile
inFile = args.inFile

figureWidth = 5
figureHeight = 2

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width = 1
panel1Height = 1

panel2Width = 0.25
panel2Height = 1

panel1 = plt.axes([0.14,0.15,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical
panel2 = plt.axes([0.076,0.15,panel2Width/figureWidth,panel2Height/figureHeight])
panel3 = plt.axes([0.14,0.685,panel2Height/figureWidth,panel2Width/figureHeight])


panel1.set_xlim(0,15)
panel1.set_ylim(0,15)

panel2.set_xlim(20,0)
panel2.set_ylim(0,15)

panel3.set_xlim(0,15)
panel3.set_ylim(0,20)


panel1.tick_params(bottom = True, labelbottom = True,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

panel2.tick_params(bottom = True, labelbottom = True,
                  left = True, labelleft = True,
                  labelright = False, labeltop = False)

panel3.tick_params(bottom = False, labelbottom = False,
                  left = True, labelleft = True,
                  labelright = False, labeltop = False)


xList = []
yList = []

fileName = 'BME163_Input_Data_1.txt' #'/Users/Sarah/Desktop/BME 163/BME163_Input_Data_1.txt'
    
for line in open(inFile):
    a = line.rstrip().split('\t')
    method1 = math.log2(int(a[1])+1)
    method2 = math.log2(int(a[2])+1)
    
    xList.append(method1)
    yList.append(method2)

panel1.plot(xList,yList,
            marker = 'o',
            markersize = 1.5,
            linewidth = 0,
            linestyle = '--',
            markeredgewidth = 0,
            markeredgecolor = 'black',
            markerfacecolor = 'black',
            alpha = 0.1 #opacity
            )

bins = np.arange(0,15,0.5)

xHisto, bins = np.histogram(xList, bins)

for index in range(0,len(xHisto),1):
    bottom = 0
    left = bins[index]
    width= bins[index+1] - left
    height = xHisto[index]

    rectangle1 = mplpatches.Rectangle([left,bottom],width,math.log2(height+1),
                                 linewidth = 0.1,
                                 facecolor = 'grey',
                                 edgecolor = 'black')
    
    panel3.add_patch(rectangle1)
    
    
bins = np.arange(0,20,0.5)

yHisto, bins = np.histogram(yList, bins)

for index in range(0,len(yHisto),1):
    
    bottom = bins[index]
    left = 0
    width = yHisto[index]
    height = bins[index+1] - bottom

    rectangle2 = mplpatches.Rectangle([left,bottom],math.log2(width+1),height,
                                 linewidth = 0.1,
                                 facecolor = 'grey',
                                 edgecolor = 'black')
    
    panel2.add_patch(rectangle2)

plt.savefig(outFile, dpi = 600)



