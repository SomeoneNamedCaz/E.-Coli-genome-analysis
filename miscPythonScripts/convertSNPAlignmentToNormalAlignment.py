import sys
from functions import *
from collections import deque
snpGenomePath = sys.argv[1]
snpIndexesPath = sys.argv[2]
print("snp genome path", snpGenomePath)
print("index path", snpIndexesPath)

refGenomePath = sys.argv[3] # .gbff file
refGenomeContigs = getContigs(refGenomePath)
if len(refGenomeContigs) > 1:
    print("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")

refSeq = refGenomeContigs[0]

outputPath = sys.argv[4]  # with file

snpIndexes = []
# load indexes
with open(snpIndexesPath) as file:

    for line in file:  # special stuff here to adjust for float method that was added
        line = float(line.strip())
        snpIndexes.append(line - 1)  # because genome files count starting at 1


# load snp genomes
alignedGenomes = {}
with open(snpIndexesPath) as file:
    fileName = "somethingBroke" # should be written to first
    for line in file:
        line = line.strip()
        if line[0] == ">":
            fileName = line[1:]
        else:
            alignedGenomes[fileName] = deque(line)



lastIndex = 0
for index in snpIndexes:
    # indexIsInsertion = False # unnecessary I htink
    # if index - int(index) != 0:
    #     indexIsInsertion = True
    for key in alignedGenomes.keys():
        if not indexIsInsertion:
            alignedGenomes[key].insert(refSeq[lastIndex:index])



