from functions import *
from findNamesOfGroups import findNamesOfGroups
"""
gets significant SNPs and maps them to the reference genome
NEEDS: path of all results
       path of a file with the indexes of the SNPS (one index on each line, in order)
       
       
       
       
       
       path command:
       time python parsingMegaCatsResults.py ./megaCATS-main/k-12/deadliness-rMsaInput.txt-rResultChisqTestFixed.txt ./megaCATS-main/k-12/allSnpsIndexes.txt K-12Ref ./refGenomes/k-12.gbff 
"""

#B2 has no nonsparse significant snps so throughs an error metadata string out og range

codonsToAminoAcids = makeCodonDictionary()

if len(sys.argv) < 7:
    print("""please give arguments: combined megaCats file, file with the snp genomes, file with the snp indexes, the suffix you want for output files,
          the reference genome file, the metadata used for megaCats, optional: whether or not you want to remove the sparce entries (true by default), the output directory (current directory by default)""")
    exit(1)


print("_________")
# args
snpStatPath = sys.argv[1]#"./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCatsPathFixed.tsv.txt"
snpGenomePath = sys.argv[2]#"./InsertAndDeleteCombinedGenomes/insertAndDeleteIndexes.txt"
snpIndexesPath = sys.argv[3]
print("stat",snpStatPath)
print("genome", snpGenomePath)
print("index", snpIndexesPath)
# add FrameShiftedToIndexPath
snpIndexesFrameShiftedPath = ".".join(snpIndexesPath.split(".")[:-1]) + "FrameShifted." + snpIndexesPath.split(".")[-1]
metaDataFilePath = sys.argv[6]
# metaDataFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/mastitisSeverityMetadata.tsv" # for severity

suffix = sys.argv[4]

annotatedRefGenomePath = sys.argv[5]#"./refGenomes/k-12.gbff"#"./AllAssemblies/refGenomeAnnotationsEdited.gb"
removeSparce = True
try:
    if sys.argv[7] == "False":
        removeSparce = False
except IndexError:
    pass

try:
    outputLocation = sys.argv[8]
except IndexError:
    outputLocation = "."

if outputLocation[-1] == "/":
    outputLocation = outputLocation[:-1]

numGenesToInclude = 10000
numSnpsToIncludeForMostSigSnps = 10_000_000


# outputs
snpsFileWithCorrectPosPath = outputLocation + "/sigSNPsByPosOnRefGenome" + suffix + ".txt" # created and then read
snpsSortedBySignificancePath =  outputLocation + "/snpsSortedBySignificanceWithGenesContainingThem" + suffix # created (extension added later)
numSnpsWithinGenesPath =  outputLocation + "/numSnpsWithinGenes"


namesOfGroups = findNamesOfGroups(metaDataFilePath, snpGenomePath, snpStatPath)#{'animal': ['chicken','cow'], 'pathogenicity':['commensal', "pathogen"], 'deadliness':['commensal', "pathogen"]}
# print(nam)
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
        currMetaDataColName = cols[-1].split("-")[0]
        positionInSnpGenome = int(cols[0])
        posInRefGenome = snpLocations[positionInSnpGenome - 1]
        if pval < significanceLevel and (not removeSparce or cols[4].strip() == "N"):
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
    print("no frame shifted index file found. Frame shifts won't be analyzed")
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
                seqWithoutSNP = gene.sequence[snp.location - (snp.location % 3) : snp.location+3-(snp.location % 3)]
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
        outFile.write("SNP pValue\tGroupEnrichedInSNP\tSNPlocation\tgeneName\tnewNuc\tsnpGroupDistribution\toldNuc\tWildTypeDistribution\tmutationType(frameShiftRecordedOnlyOnFirstNucOfIndel)\tgeneSequence\n")  # add category and move to the output function
        for snpAndGene in mostSignificantSnps[:numSnpsToIncludeForMostSigSnps]:
            snp = snpAndGene[0]
            gene = snpAndGene[1]
            #                                this is the index of the group more likely to have snp
            snpGroup = snp.allNucCounts[1-namesOfGroups[metadataCategory].index(snp.nameOfGroupMoreLikelyToHaveSNP)]
            nonSnpGroup = snp.allNucCounts[namesOfGroups[metadataCategory].index(snp.nameOfGroupMoreLikelyToHaveSNP)]
            outFile.write(
                str(snpAndGene[0].pValue) + "\t" + snp.nameOfGroupMoreLikelyToHaveSNP + "\t" + str(snpAndGene[0].location) + "\t" + snpAndGene[1].name + "\t"
                + snp.newNuc + "\t" + str(snpGroup) + "\t" + snp.oldNuc + "\t" + str(nonSnpGroup) + "\t" + str(snp.mutationType.name) + "\t" + snpAndGene[1].sequence + "\n")




import cProfile
prof = cProfile.Profile()
prof.enable()


if len(snpsFileWithCorrectPosData) == 0:
    print("no significant snps found")
    exit(10)

lastMetaDataColName = ""
indexOfLastGene = 0
lastSNPpos = 0
x = 0
for line in snpsFileWithCorrectPosData:
    # print(line)
    if x % 10000 == 0:
        print(x)
    x += 1
    line = line.strip()
    cols = line.split("\t")
    if line == "" or cols[0] == "Position" or line[0] == '"':
        continue
    pval = float(cols[2])
    currMetaDataColName = cols[-1].split("-")[0].lower()
    positionInGenome = int(cols[0])
    nucInfo = cols[5]

    groups = nucInfo.split("|")
    oldNuc, newNuc, indexOfMostSnpedGroup, _ = getSnpInfo(nucInfo)


    otherNucs = set()
    numsAndNucs = getNumsOfNucs(nucInfo)
    for numAndNuc in numsAndNucs:
        otherNucs = otherNucs.union(numAndNuc.keys())
    otherNucs.remove(oldNuc)
    otherNucs.remove(newNuc)


    # group index is the index of the higher group with the snps

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

    lastMetaDataColName = currMetaDataColName
    lastSNPpos = positionInGenome
    index = indexOfLastGene - 1
    for gene in genes:#
        index += 1

        # if in right gene
        if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
            gene.counter += 1
            #   this doesn't necessarily capture the case where group1(336_A,_75_G,_49_T)|group2(367_A,_18_G,_98_T)
            # print("currMetaDataColName\n\n", "'"+ currMetaDataColName + "'", "\n\n\n")
            snpGroup = namesOfGroups[currMetaDataColName.strip()][1-indexOfMostSnpedGroup]

            gene.snps.append(SNP(positionInGenome - gene.startPos, oldNuc, newNuc, numsAndNucs, pval, snpGroup))
            break

weights = {}
for gene in genes:
    weight = 0
    for snp in gene.snps:
        weight += math.frexp(snp.pValue)[1]
    weights[gene] = weight / len(gene.sequence)
def sortFunc(geneArg):
    return weights[geneArg]
genes.sort(key=sortFunc)
# genes.sort(key=lambda gene: 1 - gene.counter/len(gene.sequence)) # to sort by ascending proportion
print(lastMetaDataColName)
outputFunction(genes, lastMetaDataColName, weights) # "Animal"


prof.disable()
prof.print_stats(sort=1)