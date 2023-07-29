#!/usr/bin/env python
# coding: utf-8

# In[1]:
# Run time is about 20 seconds

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 


plt.style.use('BME163.mplstyle')


parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')
parser.add_argument('--inFile','-i',type = str,action = 'store',help = 'input file')

args = parser.parse_args()
outFile = args.outFile
inFile = args.inFile

figureWidth = 5
figureHeight = 3

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth = 0.75
panelHeight = 2.5

panel = plt.axes([0.1,0.1,panelWidth/figureWidth,panelHeight/figureHeight]) #horizontal, vertical
plt.xlabel('CT')
plt.ylabel('Number of genes')

panel.set_xlim(0,8)
panel.set_ylim(0,1262)


panel.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
panel.set_xticklabels(labels=['0','','6','','12','','18',''])

#########################################################################################
viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)

R1=np.linspace(viridis1[0],viridis2[0],25)
G1=np.linspace(viridis1[1],viridis2[1],25)
B1=np.linspace(viridis1[2],viridis2[2],25)

R2=np.linspace(viridis2[0],viridis3[0],25)
G2=np.linspace(viridis2[1],viridis3[1],25)
B2=np.linspace(viridis2[2],viridis3[2],25)

R3=np.linspace(viridis3[0],viridis4[0],25)
G3=np.linspace(viridis3[1],viridis4[1],25)
B3=np.linspace(viridis3[2],viridis4[2],25)

R4=np.linspace(viridis4[0],viridis5[0],26)
G4=np.linspace(viridis4[1],viridis5[1],26)
B4=np.linspace(viridis4[2],viridis5[2],26)

R=np.concatenate((R1,R2,R3,R4),axis=None)
G=np.concatenate((G1,G2,G3,G4),axis=None)
B=np.concatenate((B1,B2,B3,B4),axis=None)

########################################################################

CT = {}


with open(inFile) as f:
    next(f)
    for line in f:

        a = line.rstrip().split('\t')
        
        fpkm = [int(i) for i in a[4:12]]
        peakPhase = float(a[13])

        normalized = []

        for num in fpkm:
            num = (((num)-min(fpkm)) / (max(fpkm)-min(fpkm))) * 100
            normalized.append(num)
        if peakPhase in CT:
            CT[peakPhase].append(normalized)
        else:
            CT[peakPhase] = []
            CT[peakPhase].append(normalized)

sortedCT = dict(sorted(CT.items(), key=lambda x: x[0], reverse=True))  

plotList = []

for key in sortedCT:
    for newList in sortedCT[key]:
        plotList.append(newList)
    

for index in range(len(plotList)):
    for value in range(len(plotList[index])):
        num = int(plotList[index][value])
        rectangle1 = mplpatches.Rectangle([value,index],1,1,     # (xy, width, height)
                                             linewidth=0,
                                             facecolor=(R[num],G[num],B[num])
                    )
        panel.add_patch(rectangle1)
            



plt.savefig(outFile, dpi = 600)



