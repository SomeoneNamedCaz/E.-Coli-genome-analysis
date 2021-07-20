from functions import *
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (one index on each line, in order)
       numberOfLines in snpStatPath
       numPvaluesForEachMetadataCatagory
       
"""

codonsToAminoAcids = makeCodonDictionary()

if len(sys.argv) < 5:
    print("""please give arguments: combined megaCats file, indexes of the snpsFile, the suffix you want for output files,
          the reference genome file""")
    exit(1)

# args
snpStatPath = sys.argv[1]#"./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCatsPathFixed.tsv.txt" #TODO: send fixed stuff to doctor Erickson
snpIndexesPath = sys.argv[2]#"./InsertAndDeleteCombinedGenomes/insertAndDeleteIndexes.txt"

# add FrameShiftedToIndexPath
snpIndexesFrameShiftedPath = ".".join(snpIndexesPath.split(".")[:-2] + [snpIndexesPath.split(".")[-2]+"FrameShifted"] + [snpIndexesPath.split(".")[-1]])

suffix = sys.argv[3]

# arg less likely to change
annotatedRefGenomePath = sys.argv[4]#"./refGenomes/k-12.gbff"#"./AllAssemblies/refGenomeAnnotationsEdited.gb"
numGenesToInclude = 10000
numSnpsToIncludeForMostSigSnps = 100_000


# outputs
snpsFileWithCorrectPosPath = "./sigSNPsByPosOnRefGenome" + suffix + ".txt" # created and then read
snpsSortedBySignificancePath = "./snpsSortedBySignificanceWithGenesContainingThem" + suffix # created (extension added later)
numSnpsWithinGenesPath = "./numSnpsWithinGenes"





percentSNPsCutOffForPercentSNPs = 0.02 # for pathway analysis

snpLocations = [] # [snplocation1, ...]


with open(snpIndexesPath) as file:
    for line in file:  # special stuff here to adjust for float method that was added
        line = float(line.strip())
        line = round((line - int(line)) * 1000 + int(line))
        snpLocations.append(line)




print(len(snpLocations))
significanceLevel = 0.05/len(snpLocations)  # can change initial P-value cutoff if wanted
snpsFileWithCorrectPosData = []
with open(snpStatPath) as file, open(snpsFileWithCorrectPosPath, "w") as outFile:
    for line in file:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Position" or line[0] == '"':
            continue
        pval = float(cols[2])
        comparisonGroup = cols[-1].split("-")[0]
        positionInSnpGenome = int(cols[0])
        posInRefGenome = snpLocations[positionInSnpGenome - 1]
        if pval < significanceLevel:
            dataToWrite = str(posInRefGenome) + "\t" + "\t".join(cols[1:]) + "\n"
            outFile.write(dataToWrite)
            snpsFileWithCorrectPosData.append(dataToWrite)
#         else:
#             print("not significant")
indexesOfFrameShiftSnps = set()
try:
    with open(snpIndexesFrameShiftedPath) as frameShiftFile:
        for line in frameShiftFile:
            indexesOfFrameShiftSnps.add(int(line))
except FileNotFoundError:
    pass

t1 = time.time()
contigs = getContigs(annotatedRefGenomePath)
genes = getGenesOnContigsByPosition(annotatedRefGenomePath, contigs)
print(time.time()-t1)
def outputFunction(listOfGenes, metadataCategory, weights):
    with open(numSnpsWithinGenesPath + metadataCategory[0].upper() + metadataCategory[1:] + suffix + ".tsv", "w") as outFile:
        # header line
        outFile.write("""Name\tProduct\tweight\tPercent Of Nucleotides That Are Significant SNPs\tgene Sequence\tPercent Nucleotides That Can Be Significant NonSynonymous Mutations\tNon-synonymous SNP Indexes\tSynonymous SNP Indexes\tgeneStartPositionInGenome\tfraction of genes in pathway with high number of snps\tnum frame shifted\tframshift location\n""")
        # find out if there are other genes in the pathway
        geneNameToNumInPathway = {}
        geneNameToNumSignificantlySnpedInPathway = {}
        for gene in genes:

            if not gene.name[:3] in geneNameToNumInPathway.keys():
                geneNameToNumInPathway[gene.name[:3]] = 1
            else:
                geneNameToNumInPathway[gene.name[:3]] += 1
            if len(gene.snps)/len(gene.sequence) > percentSNPsCutOffForPercentSNPs:
                if not gene.name[:3] in geneNameToNumSignificantlySnpedInPathway.keys():
                    geneNameToNumSignificantlySnpedInPathway[gene.name[:3]] = 1
                else:
                    geneNameToNumSignificantlySnpedInPathway[gene.name[:3]] += 1
        # data
        for i in range(min(numGenesToInclude, len(listOfGenes))):
            gene = listOfGenes[i]
            numNonSyn = 0
            nonSynSnpPositions = []
            synSnpPositions = []
            frameShiftPositions = []
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

                if not indexesOfFrameShiftSnps.isdisjoint({snp.location + gene.startPos}):
                    snp.mutationType = SNP.mutationType.frameShift
                    frameShiftPositions.append(str(snp.location))
                else:
                    # I actually can't do a nuc by nuc analysis of inserts and deletes to check for nonsynonymous mutations
                    try:
                        newAA = translate(seqWithSNP, codonsToAminoAcids)
                        if newAA != translate(seqWithoutSNP, codonsToAminoAcids):
                            if newAA == "*":
                                snp.mutationType = SNP.mutationType.earlyStop
                            else:
                                snp.mutationType = SNP.mutationType.missSense
                            numNonSyn += 1
                            nonSynSnpPositions.append(str(snp.location))
                        else:
                            snp.mutationType = SNP.mutationType.silent
                            synSnpPositions.append(str(snp.location))
                    except KeyError: # i.e. if it's an indel
                        snp.mutationType = SNP.mutationType.missSense # all indels are missense or frameshift
                        pass
            # print(gene.counter, len(gene.snps))
            outFile.write(gene.name + "\t" + gene.product + "\t" + str(weights[gene]) + "\t" + str(gene.counter /
                          len(gene.sequence)) + "\t" +
                          gene.sequence + "\t" + str(numNonSyn / len(gene.sequence)) + "\t" +
                          " ".join(nonSynSnpPositions) + "\t" + " ".join(synSnpPositions) +
                          "\t" + str(gene.startPos) + "\t")
            try:
                outFile.write(str(geneNameToNumSignificantlySnpedInPathway[gene.name[:3]])
                              + "/" + str(geneNameToNumInPathway[gene.name[:3]]) + "\t")
            except KeyError: # if no sig snps in pathway
                outFile.write("0/" + str(geneNameToNumInPathway[gene.name[:3]]) + "\t")
            outFile.write(str(len(frameShiftPositions) / len(gene.sequence)) + "\t" +
                          " ".join(frameShiftPositions) + "\n")


    """ This segment looks at the most significant snps and the genes that contain them
    """
    mostSignificantSnps = []  # sorted, most significance first
    for gene in listOfGenes:
        for snp in gene.snps:
            mostSignificantSnps.append((snp, gene))
    mostSignificantSnps.sort(key=lambda snpAndGene: snpAndGene[0].pValue)
    with open(snpsSortedBySignificancePath + metadataCategory[0].upper() + metadataCategory[1:] + ".tsv", "w") as outFile:
        outFile.write("SNP pValue\tSNPlocation\tgeneName\tnewNuc\toldNuc\tmutationType(frameShiftRecordedOnlyOnFirstNucOfIndel)\tgeneSequence\n")  # add category and move to the output function
        for snpAndGene in mostSignificantSnps[:numSnpsToIncludeForMostSigSnps]:
            snp = snpAndGene[0]
            gene = snpAndGene[1]
            outFile.write(
                str(snpAndGene[0].pValue) + "\t" + str(snpAndGene[0].location) + "\t" + snpAndGene[1].name + "\t"
                + snp.newNuc + "\t" + snp.oldNuc + "\t" + str(snp.mutationType.name) + "\t" + snpAndGene[1].sequence + "\n")




import cProfile
prof = cProfile.Profile()
prof.enable()



lastMetaDataColName = ""
indexOfLastGene = 0
lastSNPpos = 0
x = 0
for line in snpsFileWithCorrectPosData:
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
            break # so we don't say the wild type is the second and first most common nuc


    if positionInGenome + 3000 < lastSNPpos: # if it went backwards by a lot then we know we are starting a new metadata category
        print("outputting first metadata category")
        print(positionInSnpGenome, lastSNPpos)
        print(lastMetaDataColName)
        print("_-----")
        indexOfLastGene = 0
        # genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence))
        weights = {}
        for gene in genes:
            weight = 0
            for snp in gene.snps:
                weight += math.frexp(snp.pValue)[1]
            weights[gene] = weight / len(gene.sequence)
        def sortFunc(geneArg):
            return weights[geneArg]
        genes.sort(key=sortFunc)
        # genes.reverse()
        outputFunction(genes, lastMetaDataColName,weights)
        # reset counts
        for gene in genes:
            gene.counter = 0
            gene.snps = []

    lastMetaDataColName = comparisonGroup
    lastSNPpos = positionInGenome
    index = indexOfLastGene-1
    for gene in genes:#
        index += 1

        # if in right gene
        if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
            gene.counter += 1
            #   this doesn't necessarily capture the case where group1(336_A,_75_G,_49_T)|group2(367_A,_18_G,_98_T)
            gene.snps.append(SNP(positionInGenome - gene.startPos, highestNuc, secondHighestNuc, pval))
            break

weights = {}
for gene in genes:
    weight = 0
    for snp in gene.snps:
        weight += math.frexp(snp.pValue)[1]
    weights[gene] = weight / len(gene.sequence)
# def sortFunc(geneArg):
#     return weights[geneArg]
# genes.sort(key=sortFunc)
genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence)) # to sort by assending proportion
outputFunction(genes, lastMetaDataColName, weights) # "Animal"



prof.disable()
prof.print_stats(sort=1)