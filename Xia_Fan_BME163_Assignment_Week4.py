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

figureWidth = 10
figureHeight = 3

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width = 6
panel1Height = 2

panel1 = plt.axes([0.1,0.15,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical

x = ([1,2,3,4,5,6,7,8,9,10,11])
plt.xticks(x, labels=[1,2,3,4,5,6,7,8,9,10,'>10'])

panel1.set_yticks([75,80,85,90,95,100])

panel1.set_xlim(0,12)
panel1.set_ylim(75,100)


plt.xlabel("Subread Coverage")
plt.ylabel("Identity (%)")


myDict = {1:[],2:[],3:[],4:[],5:[],6:[],
         7:[],8:[],9:[],10:[],11:[]}


for line in open(inFile):
    a = line.rstrip().split('\t')
    subread = a[0]
    identity = float(a[1])
    
    n = 0
    index = 0
    place1 = 0
    place2 = 0
    
    for char in subread:
        index += 1
        if char == '_':
            n += 1
            if n == 3:
                place1 = index
            if n == 4:
                place2 = index
                
    subreadCoverage = int(subread[place1:place2-1])
    
    if subreadCoverage > 10:
        if len(myDict[11]) < 1000:
            myDict[11].append(identity)
    else:
        if len(myDict[subreadCoverage]) < 1000:
            myDict[subreadCoverage].append(identity)

            
def swarmplots(xcoord, yList):
    minDist = 1/72
    increment = minDist/5
    y_points = []
    x_points = []
    placed_points = []
    direction = 0
    stop = False
    
    for ycoord in yList:  
        placed = False
        
        if stop == True: # once the 0.4 range is reached, stop the loop
            break
        
        if len(placed_points) == 0:
            x_points.append(xcoord)
            y_points.append(ycoord)
            placed_points.append((xcoord,ycoord))
            placed = True
        else:
            for shift in np.arange(0,0.4,increment):
                if placed == False:
                    if direction%2 == 0:
                        xcoordShift = xcoord + shift
                    else:
                        xcoordShift = xcoord - shift
                    distList = []
                    for point in placed_points: 
                        x_dist = ((abs(point[0]-xcoordShift))/12)*6
                        y_dist = ((abs(point[1]-ycoord))/25)*2
                        dist = np.sqrt((x_dist**2)+(y_dist**2))
                        distList.append(dist)
                    if min(distList) > minDist:
                        x_points.append(xcoordShift)
                        y_points.append(ycoord)
                        placed_points.append((xcoordShift,ycoord))
                        placed = True
                    if (0.4 - shift) < increment: # once the 0.4 range is reached, stop the loop
                        stop = True
                        break
                           
                direction += 1
    
                   
    panel1.plot(x_points,y_points,
                    marker = 'o',
                    markersize = 1,
                    linewidth = 0,
                    linestyle = '--',
                    markeredgewidth = 0,
                    markerfacecolor = 'black'
                    )    
        
    median = np.median(yList)
    x_median = [xcoord-0.4,xcoord+0.4]
    y_median = [median,median]
    
    panel1.plot(x_median,y_median,
           linestyle = '-',
           linewidth = 0.75,
           color = 'red')
    
    print(1000-len(placed_points),'points could not be plotted at position',xcoord)   

    
for number in range(1,12):
    swarmplots(number, myDict[number])
      


plt.savefig(outFile, dpi = 600)