from concurrent.futures import *
try:
    from .functions import *
    from .findNamesOfGroups import findNamesOfGroups
except ImportError:
    from functions import *
    from findNamesOfGroups import findNamesOfGroups

def loadIndexes(indexFilePath):
    snpLocations = []  # [snplocation1, ...]
    with open(indexFilePath) as file:
        for line in file:  # special stuff here to adjust for float method that was added
            line = float(line.strip())
            line = round((line - int(line)) * 1000 + int(line))
            snpLocations.append(line)
    return snpLocations
    
def calcNumGenomesWithoutGene(refGenes, pathOfAnnotatedScaffolds, metadataPath):
    numGenomesWithoutGene = {} # [gene][metadata] = numberOfgenes
    scaffoldFutures = []
    scaffoldGenesFromAllFiles = {}
    pool = ProcessPoolExecutor(8)
    for file in glob(pathOfAnnotatedScaffolds):
        
        def x(file):
            return getGenesOnContigs(file, getContigs(file))
        
        scaffoldFutures.append((file, pool.submit(x, file)))
        
    genomeNameToMetadata = readMetaDataAsDict(metadataPath)
    print("keys",genomeNameToMetadata.keys())
    for gene in refGenes:
        numGenomesWithoutGene[gene.name] = {}
        for metadataList in genomeNameToMetadata.values():
            for value in metadataList:
                numGenomesWithoutGene[gene.name][value] = 0
    
    for future in scaffoldFutures:
        genomeName = future[0].split("/")[-1]
        metadataForGenome = genomeNameToMetadata[re.sub("(.fasta)|(.gbk)","",genomeName) + ".fasta.vcf"]
        scafGenes = future[1].result()
        for gene in refGenes:
                try:
                    scafGenes[gene.name]
                except KeyError:
                    for metadata in metadataForGenome:
                        numGenomesWithoutGene[gene.name][metadata] += 1
    return numGenomesWithoutGene
def getMutationType(gene, snp,indexesOfFrameShiftSnps):
    if gene.isForward:
        seqWithoutSNP = gene.sequence[snp.location - (snp.location % 3): snp.location + 3 - (snp.location % 3)]
    else:
        seqWithoutSNP = gene.sequence[len(gene.sequence) - (snp.location - (snp.location % 3)) - 3:
                                                         len(gene.sequence) - (snp.location - (snp.location % 3))]
    # seqWithoutSNP[snp.location] = snp.oldNuc
    seqWithoutSNP = list(seqWithoutSNP)
    seqWithSNP = deepcopy(seqWithoutSNP)
    if gene.isForward:
        seqWithSNP[snp.location % 3] = snp.newNuc
    else:
        seqWithSNP[2 - snp.location % 3] = reverseComplement(snp.newNuc)
    
    seqWithSNP = "".join(seqWithSNP)  # because strings are immutable
    
    if gene.isForward:
        seqWithoutSNP[snp.location % 3] = snp.oldNuc  # in case the ref has the snp too
    else:
        seqWithoutSNP[2 - snp.location % 3] = reverseComplement(snp.oldNuc)
    seqWithoutSNP = "".join(seqWithoutSNP)  # because strings are immutable
    
    newAA = "n/a"
    oldAA = "n/a"
    if not indexesOfFrameShiftSnps.isdisjoint({snp.location + gene.startPos}):
        snp.mutationType = SNP.mutationType.frameShift
    else:
        # I actually can't do a nuc by nuc analysis of inserts and deletes to check for nonsynonymous mutations
        try:
            newAA = translate(seqWithSNP)
            oldAA = translate(seqWithoutSNP)
            if newAA != oldAA:
                if newAA == "*":
                    snp.mutationType = SNP.mutationType.earlyStop
                else:
                    snp.mutationType = SNP.mutationType.misSense
            else:
                snp.mutationType = SNP.mutationType.silent
        except KeyError:  # i.e. if it's an indel
            snp.mutationType = SNP.mutationType.misSense  # all indels are missense or frameshift
    return snp.mutationType, oldAA, newAA


def parseMegaCatsFile(megaCatsFile, snpGenomePath, snpIndexesPath, suffix, metaDataFilePath, annotatedRefGenomePath, removeSparce=True, outputDirectory="./", pathOfAnnotatedScaffolds="/Users/cazcullimore/dev/data/k-12RefGenomeAnalysisFiles/AllAssemblies/allBovineScaffolds/*.gbk", ignoreAnnotations=False, debug=False):
    if outputDirectory[-1] != "/":
        outputDirectory += "/" #TODO:fix for windows
    # add FrameShiftedToIndexPath
    snpIndexesFrameShiftedPath = ".".join(snpIndexesPath.split(".")[:-1]) + "FrameShifted." + snpIndexesPath.split(".")[-1]
    # hard-coded parameters for output
    numGenesToInclude = 10000
    numSnpsToIncludeForMostSigSnps = 10_000_000
    percentSNPsCutOffForPercentSNPs = 0.02 # for pathway analysis
    snpLocations = loadIndexes(snpIndexesPath)
    significanceLevel = 0.05/len(snpLocations)  # can change initial P-value cutoff if wanted for manhattan plot
    # outputs
    
    snpsSortedBySignificancePath = outputDirectory + "snpsSortedBySignificanceWithGenesContainingThem" + suffix # created (extension added later)
    numSnpsWithinGenesPath = outputDirectory + "numSnpsWithinGenes"
    
    def writeOutputToFiles(listOfGenes, metadataCategory, snpsForCurrentMetadataCategory, weights,numGenesToInclude,numGenomesWithoutGene, namesOfGroups):
        # make file that looks at the most snpped genes
        with open(numSnpsWithinGenesPath + metadataCategory[0].upper() + metadataCategory[1:] + suffix + ".tsv",
                  "w") as outFile:
            # header line
            outFile.write(
                """Name\tProduct\tweight\tPercent Of Nucleotides That Are Significant SNPs\tgene Sequence\tPercent Nucleotides That Can Be Significant NonSynonymous Mutations\tNon-synonymous SNP Indexes\tSynonymous SNP Indexes\tgeneStartPositionInGenome\tfraction of genes in pathway with high number of snps\tnum frame shifted\tframshift location\n""")
            # find out if there are other genes in the pathway
            geneNameToNumInPathway = {}
            geneNameToNumSignificantlySnpedInPathway = {}
            for gene in listOfGenes:
                
                if not gene.name[:3] in geneNameToNumInPathway.keys():
                    geneNameToNumInPathway[gene.name[:3]] = 1
                else:
                    geneNameToNumInPathway[gene.name[:3]] += 1
                if len(gene.snps) / len(gene.sequence) > percentSNPsCutOffForPercentSNPs:
                    if not gene.name[:3] in geneNameToNumSignificantlySnpedInPathway.keys():
                        geneNameToNumSignificantlySnpedInPathway[gene.name[:3]] = 1
                    else:
                        geneNameToNumSignificantlySnpedInPathway[gene.name[:3]] += 1
            # data
            for i in range(min(numGenesToInclude, len(listOfGenes))):
                gene = listOfGenes[i]
                nonSynSnpPositions = []
                synSnpPositions = []
                frameShiftPositions = []
                
                for snp in gene.snps:
                    snp.mutationType = getMutationType(gene, snp, indexesOfFrameShiftSnps=indexesOfFrameShiftSnps)[0]
                    if snp.mutationType == SNP.mutationType.frameShift:
                        frameShiftPositions.append(str(snp.location))
                        nonSynSnpPositions.append(str(snp.location))
                    elif snp.mutationType == SNP.mutationType.silent:
                        synSnpPositions.append(str(snp.location))
                    else:
                        nonSynSnpPositions.append(str(snp.location))
                
                outFile.write(gene.name + "\t" +
                              gene.product + "\t" +
                              str(weights[gene]) + "\t" +
                              str(gene.counter / len(gene.sequence)) + "\t" +
                              gene.sequence + "\t" +
                              str(len(nonSynSnpPositions) / len(gene.sequence)) + "\t" +
                              " ".join(nonSynSnpPositions) + "\t" +
                              " ".join(synSnpPositions) + "\t" + str(gene.startPos) + "\t")
                try:
                    outFile.write(str(geneNameToNumSignificantlySnpedInPathway[gene.name[:3]])
                                  + "/" + str(geneNameToNumInPathway[gene.name[:3]]) + "\t")
                except KeyError:  # if no sig snps in pathway
                    outFile.write("0/" + str(geneNameToNumInPathway[gene.name[:3]]) + "\t")
                outFile.write(str(len(frameShiftPositions) / len(gene.sequence)) + "\t" +
                              " ".join(frameShiftPositions) + "\n")
        
        """ This segment looks at the most significant snps and the genes that contain them
        """
        
        snpsForCurrentMetadataCategory.sort(key=lambda snp: snp.pValue)
        with open(snpsSortedBySignificancePath + metadataCategory[0].upper() + metadataCategory[1:] + ".tsv",
                  "w") as outFile:
            outFile.write(
                "SNP pValue\t"
                "GroupEnrichedInSNP\t"
                "SNPdistanceFromGeneStartPos\t"
                "UncorrectedIndexInGene\t"
                "geneName\t"
                "newNuc\t"
                "snpGroupDistribution\t"
                "numGenomesMissingGene\t"
                "oldNuc\t"
                "WildTypeDistribution\t"
                "mutationType(frameShiftRecordedOnlyOnFirstNucOfIndel)\t"
                "IndexOfGeneStart\t"
                "geneSequence\n")  # add category and move to the output function
            for snp in snpsForCurrentMetadataCategory[:numSnpsToIncludeForMostSigSnps]:
                snpIndexInGene = snp.location
                if not snp.geneContainingSnp.isForward:
                    snpIndexInGene = len(snp.geneContainingSnp.sequence) - snp.location
                #                                this is the index of the group more likely to have snp
                snpGroupDistribution = snp.allNucCounts[
                    namesOfGroups[metadataCategory].index(snp.nameOfGroupMoreLikelyToHaveSNP)]
                nonSnpGroupDistribution = snp.allNucCounts[
                    1 - namesOfGroups[metadataCategory].index(snp.nameOfGroupMoreLikelyToHaveSNP)]
                genomesWithoutGene = "no gene found"
                try:
                    genomesWithoutGene = str(numGenomesWithoutGene[snp.geneContainingSnp.name.split(":")[-1]])
                except KeyError:
                    pass
                if snp.mutationType != SNP.mutationType.silent:
                    outFile.write(
                        str(snp.pValue) + "\t" +
                        snp.nameOfGroupMoreLikelyToHaveSNP + "\t" +
                        str(snp.location) + "\t" +
                        str(snpIndexInGene) + "\t" + snp.geneContainingSnp.name + "\t" +
                        snp.newNuc + "\t" +
                        str(snpGroupDistribution) + "\t" +
                        genomesWithoutGene + "\t" +
                        snp.oldNuc + "\t" +
                        str(nonSnpGroupDistribution) + "\t" +
                        str(snp.mutationType.name) + "\t" +
                        str(snp.geneContainingSnp.startPos) + "\t" +
                        snp.geneContainingSnp.sequence + "\n")




    namesOfGroups = findNamesOfGroups(megaCatsFile)#{'animal': ['chicken','cow'], 'pathogenicity':['commensal', "pathogen"], 'deadliness':['commensal', "pathogen"]}

    snpsForCurrentMetadataCategory = []

    if abs(len(getFirstDataLine(snpGenomePath)) - len(snpLocations)) > 10_000:
        raise Exception("indexes are probably off")
    
    
    if debug:
        print("length of snp genomes", len(getFirstDataLine(snpGenomePath).strip()))
        print("length of snp indexes",len(snpLocations))

    snpsFileWithCorrectPosData = []
    with open(megaCatsFile) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or cols[0] == "Position" or line[0] == '"':
                continue
            pval = float(cols[2])
            currMetaDataColName = cols[-1].split("-")[0]
            positionInSnpGenome = int(cols[0])
            try:
                posInRefGenome = snpLocations[positionInSnpGenome - 1]
            except IndexError:
                print(positionInSnpGenome-1, "is too far")
            if pval < significanceLevel and (not removeSparce or cols[4].strip() == "N"):
                dataToWrite = str(posInRefGenome) + "\t" + "\t".join(cols[1:]) + "\n"
                snpsFileWithCorrectPosData.append(dataToWrite)
        print("last megaCats pos",positionInSnpGenome)

    indexesOfFrameShiftSnps = set()
    try:
        with open(snpIndexesFrameShiftedPath) as frameShiftFile:
            for line in frameShiftFile:
                indexesOfFrameShiftSnps.add(int(line))
    except FileNotFoundError:
        print("no frame shifted index snp file found. Frame shifts won't be analyzed")
        pass

    t1 = time.time()
    contigs = getContigs(annotatedRefGenomePath)
    genes = getGenesOnContigsByPosition(annotatedRefGenomePath, contigs)
    if not ignoreAnnotations:
        numGenomesWithoutGene = calcNumGenomesWithoutGene(genes,pathOfAnnotatedScaffolds)
    else:
        numGenomesWithoutGene = {gene.name: "unknown" for gene in genes}
    print(numGenomesWithoutGene)
    print(time.time()-t1)

    # import cProfile
    # prof = cProfile.Profile()
    # prof.enable()


    if len(snpsFileWithCorrectPosData) == 0:
        print("no significant snps found")
        raise Exception("no data found")

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

        oldNuc, newNuc, indexOfMostSnpedGroup = getSnpInfo(nucInfo)

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
            
            weights = {}
            
            # calculate weight for genes based on number of snps and snp significance
            for gene in genes:
                weight = 0
                for snp in gene.snps:
                    weight += math.frexp(snp.pValue)[1]
                weights[gene] = weight / len(gene.sequence)
            def sortFunc(geneArg):
                return weights[geneArg]
            genes.sort(key=sortFunc)
            
            writeOutputToFiles(genes,lastMetaDataColName,snpsForCurrentMetadataCategory,weights, numGenesToInclude,numGenomesWithoutGene, namesOfGroups)
            # reset counts
            for gene in genes:
                gene.counter = 0
                gene.snps = []
            snpsForCurrentMetadataCategory = []

        lastMetaDataColName = currMetaDataColName
        lastSNPpos = positionInGenome
        index = indexOfLastGene - 1
        snpFoundWithinGene = False

        minDistance = 1e30
        nearestGene = ""
        snpGroup = namesOfGroups[currMetaDataColName.strip()][indexOfMostSnpedGroup]
        for gene in genes:#
            index += 1
            # if in right gene
            if gene.stopPos > positionInGenome and gene.startPos < positionInGenome:
                gene.counter += 1
                newSnp = SNP(positionInGenome - gene.startPos, oldNuc, newNuc, numsAndNucs, pval, snpGroup, gene)
                gene.snps.append(newSnp)
                snpsForCurrentMetadataCategory.append(newSnp)
                snpFoundWithinGene = True
                break
            newDistance = min(abs(gene.startPos - positionInGenome), abs(gene.stopPos - positionInGenome))
            if newDistance < minDistance:
                nearestGene = gene
                minDistance = newDistance
        if not snpFoundWithinGene:

            nearestGeneCopy = deepcopy(nearestGene)
            nearestGeneCopy.name = "NearestGeneIs:" + nearestGeneCopy.name
            snp = SNP(positionInGenome - nearestGene.startPos, oldNuc, newNuc, numsAndNucs, pval, snpGroup,nearestGeneCopy)
            snp.mutationType = SNP.mutationType.notWithinAGene
            snpsForCurrentMetadataCategory.append(snp)

    weights = {}
    for gene in genes:
        weight = 0
        for snp in gene.snps:
            weight += math.frexp(snp.pValue)[1]
        weights[gene] = weight / len(gene.sequence)
    def sortFunc(geneArg):
        return weights[geneArg]
    genes.sort(key=sortFunc)
    print(lastMetaDataColName)
    writeOutputToFiles(genes, lastMetaDataColName, snpsForCurrentMetadataCategory, weights, numGenesToInclude, numGenomesWithoutGene,namesOfGroups)

if __name__ == "__main__":
    """
    gets significant SNPs and maps them to the reference"""

    # B2 has no nonsparse significant snps so throughs an error metadata string out og range


    if len(sys.argv) < 7:
        print("""please give arguments: combined megaCats phylogroupSnpFile, phylogroupSnpFile with the snp genomes, phylogroupSnpFile with the snp indexes, the suffix you want for output files,
              the reference genome phylogroupSnpFile, the metadata used for megaCats, optional: whether or not you want to remove the sparce entries (true by default),
               the output directory (current directory by default)""")
        exit(1)

    print("_________")
    # args
    snpStatPath = sys.argv[
        1]  # "./megaCATS-main/firstMegaCatsRunOnData/1AA_Summary_Results_ALL-metaDataForMetaCatsPathFixed.tsv.txt"
    snpGenomePath = sys.argv[2]  # "./InsertAndDeleteCombinedGenomes/insertAndDeleteIndexes.txt"
    snpIndexesPath = sys.argv[3]
    print("stat", snpStatPath)
    print("genome", snpGenomePath)
    print("index", snpIndexesPath)

    metaDataFilePath = sys.argv[6]

    suffix = sys.argv[4]

    annotatedRefGenomePath = sys.argv[5]  # "./refGenomes/k-12.gbff"#"./AllAssemblies/refGenomeAnnotationsEdited.gb"
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
    parseMegaCatsFile(snpStatPath, snpGenomePath, snpIndexesPath, suffix, metaDataFilePath, annotatedRefGenomePath,
                      removeSparce=removeSparce, outputDirectory=outputLocation)
