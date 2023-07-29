import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse 
import matplotlib.image as mpimg


plt.style.use('BME163.mplstyle')

parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o',type = str,action = 'store',help = 'output file')

parser.add_argument('--inFile1','-i1',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile2','-i2',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile3','-g',type = str,action = 'store',help = 'input file')
parser.add_argument('--inFile4','-c',type = str,action = 'store',help = 'input file')

args = parser.parse_args()
outFile = args.outFile
inFile1 = args.inFile1
inFile2 = args.inFile2
inFile3 = args.inFile3
inFile4 = args.inFile4



figureWidth = 5
figureHeight = 5

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width = 4.9
panel1Height = 1.5

panel2Width = 4.9
panel2Height = 1.5

panel3Width = 4.9
panel3Height = 1.5

panel1 = plt.axes([0.01,0.01,panel1Width/figureWidth,panel1Height/figureHeight]) #horizontal, vertical
panel2 = plt.axes([0.01,0.33,panel2Width/figureWidth,panel2Height/figureHeight])
panel3 = plt.axes([0.01,0.65,panel3Width/figureWidth,panel3Height/figureHeight])

userInput = inFile4
inputChr,start_end = userInput.split(':')
inputStart,inputEnd = [int(i) for i in start_end.split('-')]

panel1.set_xlim(inputStart-1,inputEnd+1)
panel2.set_xlim(inputStart-1,inputEnd+1)
panel3.set_xlim(inputStart-1,inputEnd+1)

panel1.tick_params(bottom = False, labelbottom = False,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

panel2.tick_params(bottom = False, labelbottom = False,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

panel3.tick_params(bottom = False, labelbottom = False,
                  left = False, labelleft = False,
                  labelright = False, labeltop = False)

##################################################################################################

def readPSL(pslFile):
    psl_list = []
    for line in open(pslFile):
        a = line.rstrip().split('\t')
        chromosome = a[13]
        start = int(a[15])
        end = int(a[16])
        blockLengths_list = [int(i) for i in a[18][:-1].split(',')]
        blockStarts_list = [int(i) for i in a[20][:-1].split(',')]
        
        if chromosome == inputChr:
            if inputStart  <= start and inputEnd >= start and inputStart  <= end and inputEnd >= end:
                psl_list.append([start,end,blockLengths_list,blockStarts_list,False])   
                
    return psl_list
                


def readGTF(gtfFile):
    gtf = {}
    gtf_list = []
    with open(inFile3) as f:
        lines_after_5 = f.readlines()[5:]
        for line in lines_after_5:
            a = line.rstrip().split('\t')
            chromosome = a[0]
            transcriptType = a[2]
            start = int(a[3])
            end = int(a[4])
            
            if chromosome == inputChr:
                if transcriptType == 'exon' or transcriptType == 'CDS':
                    if inputStart  <= start and inputEnd >= start and inputStart  <= end and inputEnd >= end:
                        b = a[8].split(';')
                        y1,transcript_id,y2 = b[1].split('"')
                        if transcript_id in gtf:
                            gtf[transcript_id].append((start,end,transcriptType))
                        else:
                            gtf[transcript_id] = [(start,end,transcriptType)]
                elif transcriptType == 'transcript':
                    if inputEnd >= start and inputStart  <= end and inputEnd >= end:
                        b = a[8].split(';')
                        y1,transcript_id,y2 = b[1].split('"')
                        if transcript_id in gtf:
                            gtf[transcript_id].append([start,end,transcriptType,False])
                        else:
                            gtf[transcript_id] = [[start,end,transcriptType,False]]
   
    alignment_list = []
    for transcript_id in gtf:
        for alignment in gtf[transcript_id]:
            alignment_list.append(alignment)
        alignment_list.sort(key = lambda x:x[0])
        gtf[transcript_id] = alignment_list
        alignment_list = []
        
    for transcript_id in gtf:
        gtf_list.append(gtf[transcript_id])
        
    gtf_list.sort(key = lambda x:x[0][1])       
    gtf_list.sort(key = lambda x:x[0][0])
    return(gtf_list)
               

def stackAlignments(alignments,panel,color1, color2, fileType, lw):
    
    if fileType == 'psl':
        yPos = 1
        yPos_list = []
        for index in range(len(alignments)): 
            previous_end = 0
            for alignment in alignments:
                if int(alignment[0]) > previous_end and alignment[4]==False:
                    plotAlignments(alignment, yPos, panel, color1, color2, fileType, lw)
                    previous_end = int(alignment[1])
                    alignment[4] = True
                    yPos_list.append(yPos)
            yPos += 1
        return max(yPos_list)
    else:
        yPos = 1
        previous_end = 0
        yPos_list = []
        for index in range(len(alignments)):
            previous_end = 0
            for alignment_list in alignments:
                if int(alignment_list[0][0]) > previous_end and alignment_list[0][3]==False:
                    plotAlignments(alignment_list, yPos, panel, color1, color2, fileType, lw)
                    previous_end = int(alignment_list[0][1])
                    alignment_list[0][3] = True
                    yPos_list.append(yPos)
            yPos += 1
        return max(yPos_list)
                  
                
def plotAlignments(alignment, yPos, panel, color1, color2, fileType, lw):
    if fileType == 'psl':
        rectangle1 = mplpatches.Rectangle([int(alignment[0]),yPos-0.025],
                                      int(alignment[1])-int(alignment[0]),0.05, 
                                      linewidth=lw, edgecolor = color2,
                                      facecolor = color1)
        panel.add_patch(rectangle1)
        
        for index in range(len(alignment[3])):
            rectangle2 = mplpatches.Rectangle([alignment[3][index],yPos-0.25], 
                                      alignment[2][index],0.5, 
                                      linewidth=lw, edgecolor = color2,
                                      facecolor = color1)
            panel.add_patch(rectangle2)
    
    if fileType == 'gtf':

        for item in alignment:
            if item[2] == 'transcript':
                rectangle1 = mplpatches.Rectangle([int(item[0]),yPos-0.025], 
                                          int(item[1])-int(item[0]),0.05, 
                                          linewidth=lw, edgecolor = color2,
                                          facecolor = color1)
                panel.add_patch(rectangle1)

            if item[2] == 'CDS':
                rectangle2 = mplpatches.Rectangle([int(item[0]),yPos-0.25], 
                                              int(item[1])-int(item[0]),0.5, 
                                              linewidth=lw, edgecolor = color2,
                                              facecolor = color1)
                panel.add_patch(rectangle2)

            if item[2] == 'exon':
                rectangle3 = mplpatches.Rectangle([int(item[0]),yPos-0.125], 
                                              int(item[1])-int(item[0]),0.25, 
                                              linewidth=lw, edgecolor = color2,
                                              facecolor = color1)
                panel.add_patch(rectangle3)
    
    
iYellow = (248/255,174/255,51/255)
iBlue=(88/255,85/255,120/255)


alignments1 = readPSL(inFile1)
alignments1.sort(key = lambda x:x[1], reverse = True)
maxY = stackAlignments(alignments1,panel1,iYellow,'grey','psl',0.02)
panel1.set_ylim(0,maxY+41)

alignments2 = readPSL(inFile2)
alignments2.sort(key = lambda x:x[1])
alignments2.sort(key = lambda x:x[0])
maxY = stackAlignments(alignments2,panel2,iBlue,'black','psl',0.05)
panel2.set_ylim(0,maxY+7.2)

alignments3 = readGTF(inFile3)
maxY = stackAlignments(alignments3,panel3,'grey','black','gtf',0.2)
panel3.set_ylim(0,maxY+1.8)

plt.savefig(outFile, dpi = 1200)
#chr7:45232000-45241000