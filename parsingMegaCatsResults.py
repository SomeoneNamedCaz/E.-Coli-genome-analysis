from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (each index on each line)
       numberOfLines in snpStatPath
       numPvaluesForEachMetadataCatagory
"""

annotatedRefGenomePath = "./AllAssemblies/refGenomeAnnotationsEdited.gb"
snpStatPath = "./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCatsPathFixed.tsv.txt" #TODO: send fixed stuff to doctor Erickson
snpLocationPath = "./substSNPcombinedGenomes/allSubstSNPsMoreThan9ActuallyWorkingMethodIndexes.txt"
snpsFileWithCorrectPosPath = "./sigSNPsByPosOnRefGenome.txt" # created and then read
snpsWithinGenesPath = "./snpsWithinGenes.txt" # created
numGenesToInclude = 1000
numPvaluesForEachMetadataCatagory = [417_115, 417_115] #TODO:not working properly (just uses first cut off and assumes that they are all the same

# import cProfile
# prof = cProfile.Profile()
# prof.enable()

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

t1 = time.time()
contigs = getContigs(annotatedRefGenomePath)
genes = getGenesOnContigsByPosition(annotatedRefGenomePath, contigs)
print(time.time()-t1)
def outputFunction(listOfGenes, outFileName):
    with open(outFileName, "w") as outFile:
        # header line
        outFile.write("""Name\tProduct\tPercent Of Nucleotides That Are Significant SNPs\tgene Sequence\tPercent Nucleotides That Can Be Significant NonSynonymous Mutations\tNon-synonymous SNP Indexes\tSynonymous SNP Indexes\tgeneStartPositionInGenome\n""")
        # data
        for i in range(min(numGenesToInclude, len(listOfGenes))):
            gene = listOfGenes[i]
            numNonSyn = 0
            nonSynSnpPositions = []
            synSnpPositions = []
            # print(gene.snps)
            for snp in gene.snps:

                # print("positionInGenome - gene.startPos", snp.location)
                # print(len(gene.sequence))
                seqWithoutSNP = gene.sequence[snp.location - (snp.location % 3) :snp.location+3-(snp.location % 3)]
                # seqWithoutSNP[snp.location] = snp.oldNuc
                seqWithoutSNP = list(seqWithoutSNP)
                seqWithSNP = seqWithoutSNP
                seqWithSNP[snp.location % 3] = snp.newNuc
                seqWithSNP = "".join(seqWithSNP)  # because strings are immutable
                seqWithoutSNP[snp.location % 3] = snp.oldNuc # in case if ref has the snp too
                seqWithoutSNP = "".join(seqWithoutSNP)  # because strings are immutable

                if translate(seqWithSNP) != translate(seqWithoutSNP):
                    numNonSyn += 1
                    nonSynSnpPositions.append(str(snp.location))
                else:
                    synSnpPositions.append(str(snp.location))
            # print(gene.counter, len(gene.snps))
            outFile.write(gene.name + "\t" + gene.product + "\t" + str(gene.counter/len(gene.sequence)) + "\t" + gene.sequence
                + "\t" + str(numNonSyn/len(gene.sequence)) + "\t" + " ".join(nonSynSnpPositions) + "\t" + " ".join(synSnpPositions) + "\t" + str(gene.startPos) + "\n")

lastMetaDataColName = ""
with open(snpsFileWithCorrectPosPath) as snpsFileWithCorrectPos:
    indexOfLastGene = 0
    lastSNPpos = 0
    x = 0
    for line in snpsFileWithCorrectPos:
        if x % 10000 == 0:
            print(x)
        x += 1
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position" or line[0] == '"':
            continue
        pval = float(cols[2])
        comparisonGroup = cols[-1].split("-")[0]
        positionInGenome = int(cols[0])
        nucInfo = cols[5]
        # nucInfo = nucInfo.strip()
        # nucInfo = nucInfo.split("|")
        highestNum = 0
        highestNuc = ""
        secondHighestNum = 0
        secondHighestNuc = ""


        for side in nucInfo.split("|"):

            nucAndNums = re.split("[(,)]", side)[1:-1] # cut off perentheses
            # print(nucAndNums)
            for nucAndNum in nucAndNums:
                nucAndNum = re.sub("_","",nucAndNum)
                nuc = nucAndNum[-1]
                num = int(nucAndNum[:-1])
                if num > highestNum:
                    secondHighestNuc = highestNuc
                    secondHighestNum = highestNum
                    highestNum = num
                    highestNuc = nuc
                elif num > secondHighestNum:
                    secondHighestNum = num
                    secondHighestNuc = nuc
            if len(nucAndNums) > 1:
                break # so we don't say the non variant is the second and first most common nuc
        # oldNuc = nucInfo[2]
        # newNuc = nucInfo[4]

        # side = ""  # longest side
        # for currSide in nucInfo.split("|"):  # get longest side
        #     if len(currSide) > len(side):  # just > gives 0.2685111989459816
        #         side = currSide
        # # problem found I think group1(1_A)|group2(336_A,_75_G,_49_T)|group3(367_A,_18_G,_98_T)
        # nucAndNums = re.split("[(,)]", side)[1:-1]  # cut off groups outside of perentheses
        # # print(nucAndNums)
        # for nucAndNum in nucAndNums:
        #     nucAndNum = re.sub("_", "", nucAndNum)
        #     nuc = nucAndNum[-1]
        #     num = int(nucAndNum[:-1])
        #     if num > highestNum:
        #         secondHighestNuc = highestNuc
        #         secondHighestNum = highestNum
        #         highestNum = num
        #         highestNuc = nuc
        #     elif num > secondHighestNum:
        #         secondHighestNum = num
        #         secondHighestNuc = nuc
        if positionInGenome < lastSNPpos:
            print("outputting first metadata category")
            indexOfLastGene = 0
            genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence))
            # genes.reverse()
            outputFunction(genes, lastMetaDataColName)
            # reset counts
            for gene in genes:
                gene.counter = 0
                gene.snps = []

        lastMetaDataColName = comparisonGroup
        lastSNPpos = positionInGenome
        index = indexOfLastGene-1
        for gene in genes:#[indexOfLastGene:]:
            index += 1

            # print(index)
            # if in right gene
            if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
                # print("added")
                gene.counter += 1
                #   this doesn't necessarily capture the case where group1(336_A,_75_G,_49_T)|group2(367_A,_18_G,_98_T)
                if contigs[0][positionInGenome] == secondHighestNuc:
                    print(contigs[0][positionInGenome], "snp here found", positionInGenome)
                elif contigs[0][positionInGenome] == highestNuc:
                    print("no snp here found", positionInGenome)
                else:
                    print("something else")
                gene.snps.append(SNP(positionInGenome - gene.startPos, highestNuc, secondHighestNuc))
                # if index >= 5: # saftey so does go to -1
                #
                #     indexOfLastGene = index - 5
                    # print(indexOfLastGene)
                break


genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence)) # to sort by assending proportion
outputFunction(genes, lastMetaDataColName)

# prof.disable()
# prof.print_stats(sort=1)