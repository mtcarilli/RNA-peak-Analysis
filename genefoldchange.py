#input two .fcount files, calculate fold change, log2 fold change, and subtraction for each gene


import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ip", "--immunoprecipitation",help="file path to IP file in .fcount format")
parser.add_argument("-i","--input",help="file path to input file in .fcount format")
parser.add_argument("-o","--output",help="file path to output file in .txt format")
parser.add_argument("-fc","--foldchange_cutoff",help="integer cutoff of fold change for file of enriched genes",type=float)
args = parser.parse_args()

IP = open(args.immunoprecipitation,"r")
inputfile = open(args.input,"r")

#extract genes and gene counts
IPlines = []
inputlines = []

for line in IP:
	IPlines.append(line.split("\t"))

for line in inputfile:
	inputlines.append(line.split("\t"))

#removes a descriptive line and column names
del inputlines[0]
del IPlines[0]
del inputlines[0]
del IPlines[0]

geneID = np.array(list(zip(*IPlines)[0]))
chro = np.array(list(zip(*IPlines)[1]))
start = np.array(list(zip(*IPlines)[2]))
end = np.array(list(zip(*IPlines)[3]))
strand = np.array(list(zip(*IPlines)[4]))
length = np.array(list(zip(*IPlines)[5]))

IPcounts = list(zip(*IPlines)[6])
inputcounts = list(zip(*inputlines)[6])

#removes the "\n"
for i in range(len(inputcounts)):
    inputcounts[i] = inputcounts[i].strip()

for i in range(len(IPcounts)):
    IPcounts[i] = IPcounts[i].strip()

#convert string list to numpy array of floats
IPcounts = np.array(IPcounts)
inputcounts = np.array(inputcounts)
IPcounts = IPcounts.astype(float)
inputcounts = inputcounts.astype(float)

#removes all genes with 3 or fewer UMI reads
keep = [inputcounts > 3]
keep2 = [IPcounts > 3]
geneID = geneID[keep and keep2]
chro = chro[keep and keep2]
start = start[keep and keep2]
end = end[keep and keep2]
strand = strand[keep and keep2]
lenght = length[keep and keep2]
inputcounts = inputcounts[keep and keep2]
IPcounts = IPcounts[keep and keep2]


#normalize IP and inputs by dividing by reads per million
normalizedIP = IPcounts/(np.sum(IPcounts)/1000000)
normalizedinputs = inputcounts/(np.sum(inputcounts)/1000000)

#calculate fold change, log2 fold change, and subtraction for normalized UMIs
foldchange = normalizedIP/normalizedinputs
log2foldchange = np.log2(foldchange)
subtraction = normalizedIP - normalizedinputs

#write them out!
output=open(args.output,"w")

output.write("GeneID" + "\t" + "Chr" + "\t" +"Start"+"\t"+"End"+"\t"+"Strand"+"\t"+"Length"+"\t"+"Norm IP Reads"+"\t"+"Norm Input reads"+ "\t"+"FoldChange"+"\t"+"log2FoldChange"+"\t"+"Subtract"+"\n")

for i in range(len(IPcounts)):
	output.write(geneID[i] + "\t")
	output.write(chro[i] + "\t")
	output.write(start[i] + "\t")
	output.write(end[i] + "\t")
	output.write(strand[i] + "\t")
	output.write(length[i] + "\t")
	output.write(str(normalizedIP[i]) + "\t")
	output.write(str(normalizedinputs[i]) + "\t")
	output.write(str(foldchange[i]) + "\t")
	output.write(str(log2foldchange[i]) + "\t")
	output.write(str(subtraction[i]) + "\n")
	

IP.close()
inputfile.close()
output.close()

#if the optional -fc is given, write file with genes enriched above given fold change
if (args.foldchange_cutoff == None):
    quit()
else:
    a = args.foldchange_cutoff
    enriched_genes= open("./{}_foldchange_{}".format(args.output,a),"w")
    enriched_genes.write("GeneID" + "\t" + "Chr" + "\t" +"Start"+"\t"+"End"+"\t"+"Strand"+"\t"+"Length"+"\t"+"Norm IP Reads"+"\t"+"Norm Input reads"+ "\t"+"FoldChange"+"\t"+"log2FoldChange"+"\t"+"Subtract"+"\n")
    
    
    for i in range(len(IPcounts)):
        if foldchange[i] >= a:
            enriched_genes.write(geneID[i] + "\t")
            enriched_genes.write(chro[i] + "\t")
            enriched_genes.write(start[i] + "\t")
            enriched_genes.write(end[i] + "\t")
            enriched_genes.write(strand[i] + "\t")
            enriched_genes.write(length[i] + "\t")
            enriched_genes.write(str(normalizedIP[i]) + "\t")
            enriched_genes.write(str(normalizedinputs[i]) + "\t")
            enriched_genes.write(str(foldchange[i]) + "\t")
            enriched_genes.write(str(log2foldchange[i]) + "\t")
            enriched_genes.write(str(subtraction[i]) + "\n")
        else:
            continue

    enriched_genes.close()

