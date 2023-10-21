from functions import *
def getGenomeNameWithNucAtSpecificPosition(pathOfFiles, geneName, position, nucToSearchFor):
    for filePath in glob(pathOfFiles):
        fileName = filePath.split("/")[-1]
        contigs = getContigs(filePath)
        try:
            qorBSeq = getGenesOnContigs(filePath, contigs)["qorB"][1]
            print(fileName, "has original Nuc:", qorBSeq[position] == nucToSearchFor, end="\t")
            print(qorBSeq)
        except KeyError:
            print(fileName, "has no gene named " + geneName)

from glob import glob
# mastitis
# PATH_OF_FILES = "/Users/cazcullimore/Documents/ericksonLabCode/reannotatedMastitisGenomes/gbks/*"

# bovinecommensal
# PATH_OF_FILES = "/Users/cazcullimore/Documents/ericksonLabCode/ProkkaOutBovineCommenSuperRun/*/*.gbk"

# # avian commensal
# PATH_OF_FILES = "/Users/cazcullimore/Documents/ericksonLabCode/natureStuff/secondCommensalAnnotation/gbks/*.gbk"

# avian pathogenicity
PATH_OF_FILES = "/Users/cazcullimore/Documents/ericksonLabCode/natureStuff/secondPathogenAnnotation/gbks/*.gbk"
geneName = "qorB"
position = 149
nucToSearchFor = "G"

if __name__ == "__main__":
    getGenomeNameWithNucAtSpecificPosition(PATH_OF_FILES, geneName, position, nucToSearchFor)
