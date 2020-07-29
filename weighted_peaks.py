#input bed file with peaks and text file found using samtools depth with chr, NT position, and depth
#outputs txt file with chr, start, end, strand, Total Reads, Max Depth, Relative Depth, and reads/NT for each peak

import numpy as np
 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p","--peaks",help="path to peak file in .bed format")
parser.add_argument("-d","--depths",help="path to depth file in .txt format")
parser.add_argument("-o","--output",help="path to file output in .txt format")
args = parser.parse_args()



peaks = open(args.peaks,"r")
depth = open(args.depths,"r")

#extract column data
peaklines = []
depthlines = []

for line in peaks:
    peaklines.append(line.split("\t"))

bedchr = list(zip(*peaklines)[0])
bedchrStart = list(zip(*peaklines)[1])
bedchrEnd = list(zip(*peaklines)[2])
bedchrStrand = list(zip(*peaklines)[5])


for line in depth:
    depthlines.append(line.split("\t"))
    
depthchr = list(zip(*depthlines)[0])  #make sure depthchr and bedchr have same names!
depthNT = list(zip(*depthlines)[1])   #position of depth 
depthValue = list(zip(*depthlines)[2]) 
depthValue = np.array([ele.replace('\n','') for ele in depthValue])
depthValue = depthValue.astype(np.float)

#lists to store data you want
totalReads = []
maxPeakHeights = []
peakLength = []
readsperNT = []

for i in range(len(bedchr)):
    indPeakD= []  #list to hold depths for nt of each individual peak
    for x in range(len(depthchr)):
        if (bedchr[i] == depthchr[x]) & (depthNT[x] >= bedchrStart[i]) & (depthNT[x] <= bedchrEnd[i]):
            indPeakD.append(depthValue[x])
        else:
            continue
    
    totalReads.append(sum(indPeakD))
    maxPeakHeights.append(max(indPeakD))
    peakLength.append(len(indPeakD))
    readsperNT.append(sum(indPeakD)/len(indPeakD))

relativePeakHeights = np.array(maxPeakHeights)/(min(maxPeakHeights))

#write to file 
output=open(args.output,"w")

output.write("chr")
output.write("\t")
output.write("Start")
output.write("\t")
output.write("End")
output.write("\t")
output.write("Strand")
output.write("\t")
output.write("Tot Reads")
output.write("\t")
output.write("Max Depth")
output.write("\t")
output.write("Relative Depth")
output.write("\t")
output.write("Length")
output.write("\t")
output.write("Reads/NT")
output.write("\n")
                 
for i in range(len(peaklines)):
    #for peaks, write columns:
    #chr, Start pos, End pos, strand, total reads, maxPeakDepth, relpeakDepth, len, reads/NT
    output.write(str(bedchr[i]))
    output.write("\t")
    output.write(str(bedchrStart[i]))
    output.write("\t")
    output.write(str(bedchrEnd[i]))
    output.write("\t")
    output.write(str(bedchrStrand[i]))
    output.write("\t")
    output.write(str(totalReads[i]))
    output.write("\t")
    output.write(str(maxPeakHeights[i]))
    output.write("\t")
    output.write(str(relativePeakHeights[i]))
    output.write("\t")
    output.write(str(peakLength[i]))
    output.write("\t")
    output.write(str(readsperNT[i]))
    output.write("\n")
    
output.close()
