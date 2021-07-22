from functions import *

def determineWhichGroupsTheGenomeIsIn(genomeSeq, PathToFileToLookAt):
    """

    :param genome: any snp genome that was used in the stats
    :param metadataForGenome: metadata values for the genome (might change)
    :param linesOfFileToLookAt: a file of the form of the megacats summary results all that doesn't have the indices corrected
    :return: index of group that the genome is in (0,1,2...)
    """
    groupIndices = []
    groupGenomeIsInForEachSnp = [] #[0] * linesOfFileToLookAt[0].split("\t")[5].count("|")
    with open(PathToFileToLookAt) as file:
        lastSnpIndex = -1
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == '' or line[0] == '"' or cols[5] == "NA":
                continue

            if len(groupGenomeIsInForEachSnp) == 0:
                groupGenomeIsInForEachSnp = [0] * (cols[5].count("|") + 1)
            snpIndex = int(cols[0])
            if snpIndex + 3000 < lastSnpIndex:
                groupIndices.append(argmax(groupGenomeIsInForEachSnp))
                groupGenomeIsInForEachSnp = []
            nucTheGenomeHas = genomeSeq[snpIndex - 1]
            oldNuc, newNuc, groupSnpIsIn = getSnpInfo(cols[5])

            if nucTheGenomeHas == newNuc: # if has snp]
                groupGenomeIsInForEachSnp[groupSnpIsIn] += 1

    return groupIndices + [argmax(groupGenomeIsInForEachSnp)]




if len(sys.argv) < 4:
    print("""please give arguments: combined megaCats file (sig snps by pos), the file with all the snp genomes, and the metadata file""")
    exit(1)

sigSNPsByPosPath = sys.argv[1]
genomePath = sys.argv[2]
metadataPath = sys.argv[3]
genomeSeq = getFirstDataLine(genomePath)
genomeFile = open(genomePath)
for line in genomeFile:
    genomeName = line[1:-1]
    break
genomeFile.close()

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
for index in range(len(groupToName)): # output not completely working but can do by hand
    groupToName[groupIndexes[index]] = genomeMetadata[index]
print(groupToName)


groupToNameMetadata1 = {1:'chicken', 0: 'cow'}
groupToNameMetadata2 = {1: 'commensal', 2:"pathogen"}