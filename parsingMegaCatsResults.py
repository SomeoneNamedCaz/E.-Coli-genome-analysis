from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (each index on each line)
       numberOfLines in snpStatPath
       numPvaluesForEachMetadataCatagory
"""

annotatedRefGenomePath = "./AllAssemblies/refGenomeAnnotationsEdited.gb"
snpStatPath = "./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCats.tsv.txt"
snpLocationPath = "./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt"
snpsFileWithCorrectPosPath = "./sigSNPsByPosOnRefGenome.txt" # created and then read
snpsWithinGenesPath = "./snpsWithinGenes.txt" # created
numGenesToInclude = 100
numPvaluesForEachMetadataCatagory = [417_115, 417_115] #TODO:not working properly (just uses first cut off and assumes that they are all the same


significanceLevel = 0.05/numPvaluesForEachMetadataCatagory[0] # can change initial P-value cutoff if wanted
snpLocations = [] # [snplocation1, ...]
with open(snpLocationPath) as file:
    for line in file:
        line = line.strip()
        if line == "":
            continue
        snpLocations.append(int(line))

print(len(snpLocations))
with open(snpStatPath) as file, open(snpsFileWithCorrectPosPath, "w") as outFile:
    for line in file:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position":
            continue
        pval = float(cols[2])
        comparisonGroup = cols[-1].split("-")[0]
        positionInSnpGenome = int(cols[0])
        # print("pos",positionInSnpGenome - 1)
        posInRefGenome = snpLocations[positionInSnpGenome - 1]
        if pval < significanceLevel:
            outFile.write(str(posInRefGenome) + "\t" + "\t".join(cols[1:]) + "\n")
#         else:
#             print("not significant")

contigs = getContigs(annotatedRefGenomePath)
genes = getGenesOnContigsByPosition(annotatedRefGenomePath, contigs)
def outputFunction(listOfGenes, outFileName):
    with open(outFileName, "w") as outFile:
        # header line
        outFile.write("""Name\tProduct\tPercentOfNucleotidesAreSignificantSNPs\tgeneSequence\t
        PercentNucleotidesThatCanBeSignificantNonSynonymousMutations\n""")
        # data
        for i in range(min(numGenesToInclude, len(listOfGenes))):
            outFile.write(listOfGenes[i].name + "\t" + listOfGenes[i].product + "\t" + str(listOfGenes[i].counter/len(listOfGenes[i].sequence)) + "\t" + listOfGenes[i].sequence + "\n")

lastMetaDataColName = ""
with open(snpsFileWithCorrectPosPath) as snpsFileWithCorrectPos:
    indexOfLastGene = 0
    lastSNPpos = 0

    for line in snpsFileWithCorrectPos:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position":
            continue
        pval = float(cols[2])
        comparisonGroup = cols[-1].split("-")[0]
        positionInGenome = int(cols[0])
        nucInfo = cols[5]
        # nucInfo = re.split("[(,_)]", nucInfo)
        # oldNuc = nucInfo[2]
        # newNuc = nucInfo[4]


        if positionInGenome < lastSNPpos:
            print("outputting first metadata category")
            genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence))
            # genes.reverse()
            outputFunction(genes, lastMetaDataColName)
            # reset counts
            for gene in genes:
                gene.counter = 0
                gene.snps = []
        lastMetaDataColName = comparisonGroup
        lastSNPpos = positionInGenome
        index = -1
        for gene in genes[indexOfLastGene:]:
            # print(gene.startPos, gene.stopPos)
            index += 1
            # if in right gene
            if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
                # print("added")
                gene.counter += 1
                # gene.snps.append(SNP(positionInGenome - gene.startPos, oldNuc, newNuc))
                indexOfLastGene = index
                break


genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence))
# genes.reverse()
outputFunction(genes, lastMetaDataColName)