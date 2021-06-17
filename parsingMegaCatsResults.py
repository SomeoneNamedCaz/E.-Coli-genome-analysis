from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (each index on each line)
       numberOfLines in snpStatPath
"""

annotatedRefGenomePath = "./AllAssemblies/refGenomeAnnotationsEdited.gb"
snpStatPath = "./megaCATS-main/1AA_Summary_Results_ALL-metaDataForMetaCats.tsv.txt"
snpLocationPath = "./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt"
snpsFileWithCorrectPosPath = "./sigSNPsByPosOnRefGenome.txt"

numPvalues = 535_413


# significanceLevel = 0.05/numPvalues # can change initial P-value cutoff if wanted
# snpLocations = [] # [snplocation, ...]
# with open(snpLocationPath) as file:
#     for line in file:
#         line = line.strip()
#         if line == "":
#             continue
#         snpLocations.append(int(line))
#
# print(len(snpLocations))
# with open(snpStatPath) as file, open(outFilePath, "w") as outFile:
#     for line in file:
#         line = line.strip()
#         cols = line.split("\t")
#         if line == "" or cols[0] == "Position":
#             continue
#         pval = float(cols[2])
#         comparisonGroup = cols[-1].split("-")[0]
#         positionInSnpGenome = int(cols[0])
#         # print("pos",positionInSnpGenome - 1)
#         posInRefGenome = snpLocations[positionInSnpGenome - 1]
#         if pval < significanceLevel:
#             outFile.write(str(posInRefGenome) + "\t" + "\t".join(cols[1:]) + "\n")
#         else:
#             print("not significant")

contigs = getContigs(annotatedRefGenomePath)
genes = getGenesOnContigsByPosition(annotatedRefGenomePath, contigs)
# geneCounts = [0] * len(genes)

onAnimal = False
ifOnNewSegment = False
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
        if positionInGenome == 1:
            onAnimal = not onAnimal
        # print(positionInGenome)

        # print("pos",positionInSnpGenome - 1)
        # posInRefGenome = snpLocations[positionInSnpGenome - 1]
        index = -1
        if not onAnimal:
            for gene in genes[indexOfLastGene:]:#range(len(genes)):
                # print(gene.startPos, gene.stopPos)
                index += 1
                # if in right gene
                if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
                    # print("added")
                    gene.counter += 1/len(gene.sequence)
                    indexOfLastGene = index
                    break
                # if gene.stopPos < positionInGenome:
                #     break


genes.sort(key=lambda gene: gene.counter)
genes.reverse()
for i in range(min(100, len(genes))):
    print(genes[i].name, genes[i].product, genes[i].counter, genes[i].counter * len(genes[i].sequence), genes[i].sequence)
