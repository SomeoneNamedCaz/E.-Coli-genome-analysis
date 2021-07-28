"""for fixing pathogenicity 3rd group error"""
from functions import *

def findOddGenomeOutFor(genomeSeq, PathToFileToLookAt):
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
            if line == '' or line[0] == '"' or cols[5] == "NA":# or float(cols[2]) > 0.05/500_000:
                continue
            # print(float(cols[1]))
            if len(groupGenomeIsInForEachSnp) == 0:
                groupGenomeIsInForEachSnp = [0] * (cols[5].count("|") + 1)
            snpIndex = int(cols[0])
            if snpIndex + 3000 < lastSnpIndex:
                groupIndices.append(argmax(groupGenomeIsInForEachSnp))
                groupGenomeIsInForEachSnp = []
            nucTheGenomeHas = genomeSeq[snpIndex - 1]
            numsOfNucs = getNumsOfNucs(cols[5])
            if nucTheGenomeHas in numsOfNucs[0].keys() and not nucTheGenomeHas in numsOfNucs[1].keys() and not nucTheGenomeHas in numsOfNucs[2].keys():
                return True

    return False



if len(sys.argv) < 4:
    print("""please give arguments: combined megaCats file, the file with all the snp genomes, and the metadata file""")
    exit(1)

sigSNPsByPosPath = sys.argv[1]
genomePath = sys.argv[2]
metadataPath = sys.argv[3]
for index in range(0,944):
    indexOfGenomeToGet = index # 400 is cow
    genomeSeq = open(genomePath).readlines()[1+indexOfGenomeToGet*2]#getFirstDataLine(genomePath)
    genomeName = open(genomePath).readlines()[indexOfGenomeToGet*2][1:-1]
    # genomeFile.close()

    isOddOneOut = findOddGenomeOutFor(genomeSeq, sigSNPsByPosPath)
    if isOddOneOut:
        print(genomeName) # found scaffold_1554_SS_317.fasta.vcf # last one. the error is because no newline at end of file
