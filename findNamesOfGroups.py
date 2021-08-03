from functions import *

def determineWhichGroupsTheGenomeIsIn(genomeSeq, PathToFileToLookAt):
    """

    :param genome: any snp genome that was used in the stats
    :param metadataForGenome: metadata values for the genome (might change)
    :param linesOfFileToLookAt: a file of the form of the megacats summary results all that doesn't have the indices corrected
    :return: index of group that the genome is in (0,1,2...)
    """
    groupIndices = []
    groupProbability = [] # prob part of group 0, prob part of group 1, ...
    groupGenomeIsInForEachSnp = [] #[0] * linesOfFileToLookAt[0].split("\t")[5].count("|")
    with open(PathToFileToLookAt) as file:
        lastSnpIndex = -1
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == '' or line[0] == '"' or cols[1] == "NA" or cols[0] == "Position":# or float(cols[2]) > 0.05/500_000:
                continue
            # print(float(cols[1]))
            if len(groupGenomeIsInForEachSnp) == 0:
                groupGenomeIsInForEachSnp = [0] * (cols[5].count("|") + 1)
            snpIndex = int(cols[0])
            if snpIndex + 3000 < lastSnpIndex:
                groupIndices.append(argmax(groupGenomeIsInForEachSnp))
                groupGenomeIsInForEachSnp = []
            nucTheGenomeHas = genomeSeq[snpIndex - 1]
            # print(cols[0])
            oldNuc, newNuc, groupSnpIsIn, proportion = getSnpInfo(cols[5])
            # if proportion > 0.8:
            #     continue
            # print(line)
            numsOfNucs = getNumsOfNucs(cols[5])
            for index in range(len(groupGenomeIsInForEachSnp) - 2, len(groupGenomeIsInForEachSnp)):
                # print("i",index)
                otherGroupIndex = len(groupGenomeIsInForEachSnp) - 1 -index
                try:
                    numsOfNucs[index][nucTheGenomeHas] # groupGenomeIsInForEachSnp[index] +=
                except KeyError:
                    pass
                    print(numsOfNucs, nucTheGenomeHas)
                    print("other",otherGroupIndex)
                    groupGenomeIsInForEachSnp[otherGroupIndex] += 1e6
            # if nucTheGenomeHas == oldNuc: # if doesn't has snp
            #     pass#groupGenomeIsInForEachSnp[groupSnpIsIn] += (1 - proportion)
            # else:

    print("all group counts",groupGenomeIsInForEachSnp)
    return groupIndices + [argmax(groupGenomeIsInForEachSnp)]




if len(sys.argv) < 4:
    print("""please give arguments: combined megaCats file, the file with all the snp genomes, and the metadata file""")
    exit(1)

sigSNPsByPosPath = sys.argv[1]
genomePath = sys.argv[2]
metadataPath = sys.argv[3]
indexOfGenomeToGet = 100 # 400 is cow
genomeSeq = open(genomePath).readlines()[1+indexOfGenomeToGet*2]#getFirstDataLine(genomePath)
genomeName = open(genomePath).readlines()[indexOfGenomeToGet*2][1:-1]
# genomeFile.close()

groupIndexes = determineWhichGroupsTheGenomeIsIn(genomeSeq, sigSNPsByPosPath)
genomeMetadata = []
with open(metadataPath) as metadataFile:
    for line in metadataFile:
        cols = line.strip().split("\t")
        if cols[0] == genomeName:
            genomeMetadata = cols[1:]
            break

groupToName = {}
print(groupIndexes)
print(genomeMetadata)
for index in range(len(groupIndexes)): # output not completely working but can do by hand
    groupToName[groupIndexes[index]] = genomeMetadata[index]
# print(groupToName)


groupToNameMetadata1 = {0:'chicken', 1: 'cow'}
groupToNameMetadata2 = {1: 'commensal', 2:"pathogen"}

#NOTE: first genoem is in group 2 (index 1) for animal
# or maybe not I'm not sure

# need to balance proportion by which snp it has

