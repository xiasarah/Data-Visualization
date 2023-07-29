#!/usr/bin/env python
# coding: utf-8

# In[81]:


import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 
import matplotlib.image as mpimg



plt.style.use('BME163.mplstyle')

parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')

parser.add_argument('--inFile1','-b',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile2','-g',type = str,action = 'store',help = 'input file')

parser.add_argument('--inFile3','-A',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile4','-T',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile5','-C',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile6','-G',type = str,action = 'store',help = 'input file')



args = parser.parse_args()
outFile = args.outFile
inFile1 = args.inFile1
inFile2 = args.inFile2
inFile3 = args.inFile3
inFile4 = args.inFile4
inFile5 = args.inFile5
inFile6 = args.inFile6






figureWidth = 5
figureHeight = 2

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width = 1.5
panel1Height = 0.5

panel2Width = 1.5
panel2Height = 0.5

panel1 = plt.axes([0.1,0.3,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical
plt.xlabel('Distance to\nSplice Site')
plt.ylabel('Bits')
plt.title("5'SS")
panel2 = plt.axes([0.44,0.3,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical
plt.xlabel('Distance to\nSplice Site')
plt.title("3'SS")


panel1.set_yticks([0,1,2])
panel1.set_xticks([-10,-5,0,5,10])
panel2.set_xticks([-10,-5,0,5,10])

panel1.set_xlim(-10,10)
panel2.set_xlim(-10,10)
panel1.set_ylim(0,2)
panel2.set_ylim(0,2)


panel2.tick_params(bottom = True, labelbottom = True,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

x = [0,0]
y = [0,2]
panel1.plot(x,y,
           linestyle = '-',
           linewidth = 0.5,
           color = 'black')
panel2.plot(x,y,
           linestyle = '-',
           linewidth = 0.5,
           color = 'black')

############# Fasta parser #############

seqDict = {}
name = ''
sequence = ''

for line in open(inFile2):
    if line[0] == '>':
        if name:
            seqDict[name] = ('').join(sequence)
        if line[1:].strip().split()[0] == 'MT':
            name = 'chrM'
        else:
            name = 'chr'+ line[1:].strip().split()[0]
        sequence = []
    else:
        sequence.append(line.strip())
        
if sequence:
    seqDict[name] =  ('').join(sequence)


##########################################

def revComp(seq):
    revComp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    revSeq = seq[::-1]
    seq = ''
    for nuc in range(len(revSeq)):
        seq += revComp[revSeq[nuc]]
    return seq

##########################################

fastDict = {}


nucDict5={}
nucDict3 = {}

for number in range(0,20):
    nucDict5[number] = []
    nucDict3[number] = []
    
# reading through bed file

for line in open(inFile1):
    
    a = line.rstrip().split('\t')
    chromosome = a[0]
    position = int(a[1])
    direction = a[3][:2]
    
    if direction == "5'":
        spliceSite = seqDict[chromosome][position-10:position+10]
        for nuc in range(len(spliceSite)):
            nucDict5[nuc].append(spliceSite[nuc])
            
    if direction == "3'":
        spliceSite = seqDict[chromosome][position-10:position+10]
        rc_spliceSite = revComp(spliceSite)
        for nuc in range(len(rc_spliceSite)):
            nucDict3[nuc].append(rc_spliceSite[nuc])

nucFreq5 = {}
nucFreq3 = {}

for number in range(0,20):
    nucFreq5[number] = {'A':0, 'T':0, 'C':0, 'G':0}
    nucFreq3[number] = {'A':0, 'T':0, 'C':0, 'G':0}
    

for pos in nucDict5:
    A = nucDict5[pos].count('A')/(len(nucDict5[pos]))
    T = nucDict5[pos].count('T')/(len(nucDict5[pos]))
    C = nucDict5[pos].count('C')/(len(nucDict5[pos]))
    G = nucDict5[pos].count('G')/(len(nucDict5[pos]))
    nucFreq5[pos]['A'] = A
    nucFreq5[pos]['T'] = T
    nucFreq5[pos]['C'] = C
    nucFreq5[pos]['G'] = G
    
for pos in nucDict3:
    A = nucDict3[pos].count('A')/(len(nucDict3[pos]))
    T = nucDict3[pos].count('T')/(len(nucDict3[pos]))
    C = nucDict3[pos].count('C')/(len(nucDict3[pos]))
    G = nucDict3[pos].count('G')/(len(nucDict3[pos]))
    nucFreq3[pos]['A'] = A
    nucFreq3[pos]['T'] = T
    nucFreq3[pos]['C'] = C
    nucFreq3[pos]['G'] = G
    
###########################################

seqLen5 = len(nucDict5[0])
seqLen3 = len(nucDict3[0])

en5 = (1/np.log(2))*(3/(2*seqLen5))
en3 = (1/np.log(2))*(3/(2*seqLen3))

A = mpimg.imread(inFile3)
T = mpimg.imread(inFile4)
C = mpimg.imread(inFile5)
G = mpimg.imread(inFile6)

height_A = 0
height_T = 0
height_C = 0
height_G = 0

height = []

for pos in nucFreq5:
    
    H = -(nucFreq5[pos]['A']*np.log2(nucFreq5[pos]['A']) +
        nucFreq5[pos]['T']*np.log2(nucFreq5[pos]['T']) +
        nucFreq5[pos]['C']*np.log2(nucFreq5[pos]['C']) +
        nucFreq5[pos]['G']*np.log2(nucFreq5[pos]['G']) )
        
    height_A = nucFreq5[pos]['A'] * (np.log2(4)-(H)+en5)
    height_T = nucFreq5[pos]['T'] * (np.log2(4)-(H)+en5)
    height_C = nucFreq5[pos]['C'] * (np.log2(4)-(H)+en5)
    height_G = nucFreq5[pos]['G'] * (np.log2(4)-(H)+en5)
    
    
    height = [(A,height_A), (T,height_T), (C,height_C), (G,height_G)]
    height.sort(key = lambda x:x[1])
    

    panel1.imshow(height[0][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,0,height[0][1]]) # extent = [left,right,bottom,top]
    panel1.imshow(height[1][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[0][1],height[1][1]])
    panel1.imshow(height[2][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[1][1],height[2][1]])
    panel1.imshow(height[3][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[2][1],height[3][1]])

for pos in nucFreq3:
    
    H = -(nucFreq3[pos]['A']*np.log2(nucFreq3[pos]['A']) +
        nucFreq3[pos]['T']*np.log2(nucFreq3[pos]['T']) +
        nucFreq3[pos]['C']*np.log2(nucFreq3[pos]['C']) +
        nucFreq3[pos]['G']*np.log2(nucFreq3[pos]['G']) )
        
    height_A = nucFreq3[pos]['A'] * (np.log2(4)-(H)+en3)
    height_T = nucFreq3[pos]['T'] * (np.log2(4)-(H)+en3)
    height_C = nucFreq3[pos]['C'] * (np.log2(4)-(H)+en3)
    height_G = nucFreq3[pos]['G'] * (np.log2(4)-(H)+en3)

    
    height = [(A,height_A), (T,height_T), (C,height_C), (G,height_G)]
    height.sort(key = lambda x:x[1])

    panel2.imshow(height[0][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,0,height[0][1]]) 
    panel2.imshow(height[1][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[0][1],height[1][1]])
    panel2.imshow(height[2][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[1][1],height[2][1]])
    panel2.imshow(height[3][0], aspect='auto', origin='upper', extent = [pos-10,pos-9,height[2][1],height[3][1]])

###########################################

            

plt.savefig(outFile, dpi = 600)




