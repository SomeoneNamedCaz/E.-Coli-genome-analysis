from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (each index on each line)
       numberOfLines in snpStatPath
"""

annotatedRefGenomePath = "./AllAssemblies/refGenomeAnnotationsEdited.gb"
snpStatPath = "./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCats.tsv.txt"
snpLocationPath = "./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt"
snpsFileWithCorrectPosPath = "./sigSNPsByPosOnRefGenome.txt" # created and then read
snpsWithinGenesPath = "./snpsWithinGenes.txt" # created
numGenesToInclude = 100

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
        PercentNucleotidesThatCanBeSignificantNonSynonymousMutations""")
        # data
        for i in range(min(numGenesToInclude, len(listOfGenes))):
            outFile.write(listOfGenes[i].name + "\t" + listOfGenes[i].product + "\t" + str(listOfGenes[i].counter) + "\t" + listOfGenes[i].sequence)

lastSNPpos = 0
with open(snpsFileWithCorrectPosPath) as snpsFileWithCorrectPos:
    indexOfLastGene = 0

    for line in snpsFileWithCorrectPos:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position":
            continue
        pval = float(cols[2])
        comparisonGroup = cols[-1].split("-")[0]
        positionInGenome = int(cols[0])

        if positionInGenome < lastSNPpos:
            genes.sort(key=lambda gene: gene.counter)
            genes.reverse()
            outputFunction(genes)
            # reset counts
            for gene in genes:
                gene.counter = 0

        index = -1
        for gene in genes[indexOfLastGene:]:
            # print(gene.startPos, gene.stopPos)
            index += 1
            # if in right gene
            if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
                # print("added")
                gene.counter += 1/len(gene.sequence)
                indexOfLastGene = index
                break


genes.sort(key=lambda gene: gene.counter)
genes.reverse()
outputFunction(genes)