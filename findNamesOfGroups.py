from functions import *

def determineGroupIndexes(genomeSeq, PathToFileToLookAt):
    """

    :param genome: any snp genome that was used in the stats
    :param metadataForGenome: metadata values for the genome (might change)
    :param linesOfFileToLookAt: a file of the form of the megacats summary results all that doesn't have the indices corrected
    :return: index of group that the genome is in (0,1,2...)
    """
    groupIndices = []
    groupProbability = [] # prob part of group 0, prob part of group 1, ...
    groupGenomeIsInForEachSnp = [] #[0] * linesOfFileToLookAt[0].split("\t")[5].count("|")
    metaDataCategoriesFound = []
    with open(PathToFileToLookAt) as file:
        lastSnpIndex = -1
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == '' or line[0] == '"' or cols[1] == "NA" or cols[0] == "Position":# or float(cols[2]) > 0.05/500_000:
                continue
            # print(float(cols[1]))
            metaDataCategory = cols[-1].split("-")[0]
            if not metaDataCategory in metaDataCategoriesFound:
                metaDataCategoriesFound.append(metaDataCategory)
            if len(groupGenomeIsInForEachSnp) == 0:
                groupGenomeIsInForEachSnp = [0] * (cols[5].count("|") + 1)
            snpIndex = int(cols[0])
            if snpIndex + 3000 < lastSnpIndex:
                print("all group counts", groupGenomeIsInForEachSnp)
                numGroupsWithCounts = 0
                for groupCounts in groupGenomeIsInForEachSnp:
                    if groupCounts > 0.1:
                        numGroupsWithCounts += 1
                if numGroupsWithCounts > 1:
                    raise Exception("incorrect genome or indexes being analyzed")
                groupIndices.append(argmax(groupGenomeIsInForEachSnp))
                groupGenomeIsInForEachSnp = []
            nucTheGenomeHas = genomeSeq[snpIndex - 1]
            # print(cols[0])
            oldNuc, newNuc, groupSnpIsIn, proportion = getSnpInfo(cols[5])
            # if proportion > 0.8:
            #     continue
            # print(line)
            lastSnpIndex = snpIndex
            numsOfNucs = getNumsOfNucs(cols[5])
            for index in range(0, 2):
                # print("i",index)
                otherGroupIndex = 1 - index
                try:
                    numsOfNucs[index][nucTheGenomeHas] # groupGenomeIsInForEachSnp[index] +=
                except KeyError:
                    # print(numsOfNucs, nucTheGenomeHas)
                    # print("other",otherGroupIndex)
                    groupGenomeIsInForEachSnp[otherGroupIndex] += 1e6
            # if nucTheGenomeHas == oldNuc: # if doesn't has snp
            #     pass#groupGenomeIsInForEachSnp[groupSnpIsIn] += (1 - proportion)
            # else:

    print("all group counts",groupGenomeIsInForEachSnp)
    numGroupsWithCounts = 0
    for groupCounts in groupGenomeIsInForEachSnp:
        if groupCounts > 0.1:
            numGroupsWithCounts += 1
    if numGroupsWithCounts > 1:
        raise Exception("incorrect genome or indexes being analyzed")
    return (groupIndices + [argmax(groupGenomeIsInForEachSnp)], metaDataCategoriesFound)



# genomeFile.close()
def findNamesOfGroups(metadataPath, genomePath, snpStatPath):
    indexOfGenomeToGet = 0  # 400 is cow
    try:
        genomeSeq = open(genomePath).readlines()[1 + indexOfGenomeToGet * 2]  # getFirstDataLine(genomePath)
    except IndexError:
        print(genomePath + " file not found or has fewer than 2 lines")
        exit(1)
    genomeName = open(genomePath).readlines()[indexOfGenomeToGet * 2][1:-1]
    groupIndexes, metadataCategoriesFound = determineGroupIndexes(genomeSeq, snpStatPath)
    allMetadataCategories = []
    possibleMetaDataValues = {} # {category:{possibilities}}
    genomeMetadata = []
    with open(metadataPath) as metadataFile:
        onFirstLine = True
        for line in metadataFile:
            line = line.strip()
            cols = line.split("\t")
            if line == "":
                continue
            if cols[0] == genomeName:
                genomeMetadata = cols[1:]
            if onFirstLine:
                allMetadataCategories = cols[1:] # assume name is first col
                for i in range(len(allMetadataCategories)):
                    allMetadataCategories[i] = allMetadataCategories[i].lower()
                onFirstLine = False
            else:
                i = 0
                for category in allMetadataCategories:
                    try:
                        possibleMetaDataValues[category].add(cols[i + 1])
                    except:
                        possibleMetaDataValues[category] = {cols[i + 1]}
                    i += 1

    groupToName = {}
    print("groupIndexes",groupIndexes)
    print(genomeMetadata)
    if len(metadataCategoriesFound) < len(allMetadataCategories):
        i = 0
        for category in allMetadataCategories:
            if not category in metadataCategoriesFound:
                groupIndexes.insert(i,0)
                i += 1
            i += 1



    for i in range(len(groupIndexes)): # output not completely working but can do by hand
        print(i)
        print(genomeMetadata)
        indexOfTheGroupThatTheGenomeIsIn = groupIndexes[i]
        metaDataForTheGroupTheGenomeIsIn = genomeMetadata[i]
        category = allMetadataCategories[i]
        print(metaDataForTheGroupTheGenomeIsIn)
        print(possibleMetaDataValues[category])
        possibleMetaDataValues[category].remove(metaDataForTheGroupTheGenomeIsIn)
        otherPossibleMetaDataValue = possibleMetaDataValues[category].pop() # add other element
        if indexOfTheGroupThatTheGenomeIsIn == 0:
            groupToName[category] = [metaDataForTheGroupTheGenomeIsIn, otherPossibleMetaDataValue]
        else: # if == 1
            groupToName[category] = [otherPossibleMetaDataValue, metaDataForTheGroupTheGenomeIsIn]
            # [groupIndexes[index]
    return groupToName


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            """please give arguments: combined megaCats file, the file with all the snp genomes, and the metadata file""")
        exit(1)

    snpStatPath = sys.argv[1]
    genomePath = sys.argv[2]
    metadataPath = sys.argv[3]
    print(findNamesOfGroups(metadataPath, genomePath, snpStatPath))