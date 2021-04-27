import os
from glob import glob
import re
from functions import *

prevalenceFileName = "genePrevalances.tsv" # structured as geneName + \t + count

# get most common gene number to use as proxy for number of genomes
geneCounts = []
with open(prevalenceFileName) as prevalenceFile:
    for line in prevalenceFile:
        line = line.strip()
        cols = line.split("\t")
        geneName = cols[0]
        geneCount = int(cols[1])
        geneCounts.append(geneCount)
geneCounts.sort()
assumedNumberOfStrains = geneCounts[-100]
# read through file to get percents
genes99_100 = []
genes95_99 = []
genes15_95 = []
genes0_15 = []
totalGenes = 0
with open(prevalenceFileName) as prevalenceFile:
    for line in prevalenceFile:
        line = line.strip()
        cols = line.split("\t")
        geneName = cols[0]
        geneCount = int(cols[1])
        if geneCount/assumedNumberOfStrains >= 0.99:
            genes99_100.append(geneCount / assumedNumberOfStrains)
        elif geneCount/assumedNumberOfStrains >= 0.95:
            genes95_99.append(geneCount / assumedNumberOfStrains)
        elif geneCount/assumedNumberOfStrains > 0.15:
            genes15_95.append(geneCount / assumedNumberOfStrains)
        else:
            genes0_15.append(geneCount / assumedNumberOfStrains)
        totalGenes += 1
        # print(geneName,geneCount/maxCount)
print(len(genes99_100))
print(len(genes95_99))
print(len(genes15_95))
print(len(genes0_15))
print("total genes", totalGenes)
print("assumed num strains", assumedNumberOfStrains)
