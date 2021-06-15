from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (each index on each line)
       numberOfLines in snpStatPath
"""

snpStatPath = "./megaCATS-main/1AA_Summary_Results_ALL-metaDataForMetaCats.tsv.txt"
snpLocationPath = "./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt"
outFilePath = "./sigSNPsByPosOnRefGenome.txt"

numPvalues = 535_413


significanceLevel = 0.05/numPvalues # can change initial P-value cutoff if wanted
snpLocations = [] # [snplocation, ...]
with open(snpLocationPath) as file:
    for line in file:
        line = line.strip()
        if line == "":
            continue
        snpLocations.append(int(line))

print(len(snpLocations))
with open(snpStatPath) as file, open(outFilePath, "w") as outFile:
    for line in file:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position":
            continue
        pval = float(cols[3])
        comparisonGroup = cols[-1].split("-")[0]
        positionInSnpGenome = int(cols[0])
        posInRefGenome = snpLocations[positionInSnpGenome]
        if pval < significanceLevel:
            outFile.write(str(posInRefGenome) + "\t" + "\t".join(cols[1:]) + "\n")
        else:
            print("not significant")
