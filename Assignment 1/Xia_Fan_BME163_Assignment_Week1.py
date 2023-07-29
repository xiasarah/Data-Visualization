#!/usr/bin/env python
# coding: utf-8

# In[192]:


import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 


parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')
args = parser.parse_args()
outFile = args.outFile


plt.style.use('BME163.mplstyle')

figureWidth = 3.42
figureHeight = 2

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth = 1
panelHeight = 1

panel1 = plt.axes([0.1,0.2,panelWidth/figureWidth,panelHeight/figureHeight])
panel2 = plt.axes([0.55,0.2,panelWidth/figureWidth,panelHeight/figureHeight])


# left figure

red = (1,0,0)
blue = (0,0,1)
white = (1,1,1)
black = (0,0,0)

color1 = white
color2 = black

R = np.linspace(color1[0], color2[0], 101)
G = np.linspace(color1[1], color2[1], 101)
B = np.linspace(color1[2], color2[2], 101)

for r in np.linspace(0,np.pi/2,25):

    xvalue = np.cos(r) * 100
    yvalue = np.sin(r) * 100
    
    panel1.plot(xvalue,yvalue, 
           marker = 'o',
           markersize = 2,
           linewidth = 0.5,
           linestyle = '--',
           markeredgewidth = 0,
           markerfacecolor = (R[int(xvalue)],B[int(xvalue)],G[int(xvalue)])
           )

panel1.set_xlim(0,100)
panel1.set_ylim(0,100)


# right figure

panel2.set_xlim(0,100)
panel1.set_ylim(0,100)

white = (1,1,1)
black = (0,0,0)
blue = (0,0,1)

color1 = white
color2 = blue

R = np.linspace(color1[0], color2[0], 101)
G = np.linspace(color1[1], color2[1], 101)
B = np.linspace(color1[2], color2[2], 101)

for h in np.arange(0,100,10):
    for i in np.arange(0,100,10):
        rectangle1 = mplpatches.Rectangle([i,h],10,10,  # Rectangle([left,bottom],width,height
                                 linewidth = 0.75,
                                 edgecolor='black',
                                 facecolor = (R[i],G[h],B[h]))
        panel2.add_patch(rectangle1)

    
    
panel2.set_xlim(0,100)
panel1.set_ylim(0,100)

# panel formatting

panel1.set_xlim(0,100)
panel1.set_ylim(0,100)

panel2.set_xlim(0,100)
panel2.set_ylim(0,100)
                
panel1.tick_params(bottom = False, labelbottom = False,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

panel2.tick_params(bottom = False, labelbottom = False,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)




plt.savefig(outFile, dpi = 600)


# In[ ]:




