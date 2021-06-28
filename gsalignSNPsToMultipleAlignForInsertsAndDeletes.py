"""
NOTE: this isn't working yet
"""
from glob import glob
import time
import copy
import os
import cProfile
prof = cProfile.Profile()
prof.enable()
# from functions import *
# def profileFunction():
    # gsAlignPaths = "./mastitisMultipleAlignmentStuff/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.vcf"
gsAlignPaths = "./allGsAlignOutputs/*.vcf"
# NamesOfMastitisFiles = []
snps = [] # list of elements like this: [fileName, location, oldNuc, NewNuc, type]
# snpPositions = set()
# print(snpPositions)
numFiles = 0
for filePath in glob(gsAlignPaths):
    numFiles += 1
    # print("opened one")
    with open(filePath) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or line[0] == "#":
                continue
            snps.append([filePath.split("/")[-1], int(cols[1]), cols[3], cols[4], cols[7][5:]])
            # if numFiles == 1:
            #     print([filePath, cols[1], cols[3], cols[4], cols[7][5:]])
            # snpPositions.add(cols[1])

# for snp in snps:
#     print(snp)
# print(len(snpPositions))
snps.sort(key=lambda a: a[1])
# with open("./SNPsList.tsv", "w") as file:
#     for snp in snps:
#         for val in snp:
#             file.write(val + "\t")
#         file.write("\n")




# del snpPositions
allFiles = glob(gsAlignPaths)

for index in range(len(allFiles)):
    allFiles[index] = allFiles[index].split("/")[-1]
filesWithoutSNP = copy.deepcopy(allFiles)
removedFiles = []
# print(filesWithoutSNP)
currPos = snps[0][1]
# currFile = open(outFilePathPrefix + snps[0][1] + ".afa", "w")
oldNucAtCurrentPos = snps[0][2]
# longestNucLength = determineIfVariantIncludesNumOfNucleotideChange(currPos, 0) # -1 is a flag for no insertions or deletions
snpIndex = 0
dataToWrite = {} # {fileName:SNPlist}
print(len(snps))
snpsLength = len(snps) #* len(filesWithoutSNP)
for file in allFiles: # load the list with dicts with empty strings
    dataToWrite[file] = ""
t1 = time.time()
numerrors = 0
numSkipped = 0

lastIndex = 0


maxLenOfSnpGenome = -1
snpIndex = -1
needToSkip = False
print(snps[:100])
outFilePath = "./InsertAndDeleteCombinedGenomes/" # folder
if not os.path.exists(outFilePath):
    os.mkdir(outFilePath)

# add file name
outFilePath += "insertAndDelete.afa" # needs three letter extension
pathForSNPsIncludedIndexes = outFilePath[:-4] + "Indexes.txt"

print("start loop")
with open(pathForSNPsIncludedIndexes, "w") as indexFile:
    for snp in snps:# list of elements like this: [fileName, location, oldNuc, NewNuc, type]
        snpIndex += 1
        if snp[4] != "INSERT" and snp[4] != "DELETE":
            continue
        snpPos = snp[1]

        if currPos < snpPos: # if ran out of SNPs at the position
            # add the oldNuc for that SNP
            if not needToSkip:
                for fileWithMissingSNP in filesWithoutSNP:
                    dataToWrite[fileWithMissingSNP] += oldNucAtCurrentPos
                    maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[fileWithMissingSNP]))
                for file in allFiles:
                    dataToWrite[file] += "-" * (maxLenOfSnpGenome - len(dataToWrite[file]))
                indexFile.write(str(currPos) + "\n")

            try:
                if snpPos < snps[snpIndex + 9][1]: # if less than 10 snp at current position
                    needToSkip = True
                    continue
            except:
                pass
            # reset everything for next snp
            currPos = snpPos
            filesWithoutSNP += removedFiles  # add back removed files
            removedFiles = []
            oldNucAtCurrentPos = snp[2]
        needToSkip = False
        try:
            filesWithoutSNP.remove(snp[0]) # throws error if trying to add duplicate snp
            removedFiles.append(snp[0])
            if snp[2] > len(oldNucAtCurrentPos):
                oldNucAtCurrentPos = snp[2]
            maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]]) + len(snp[2]))
            dataToWrite[snp[0]] += snp[3] # add snp to entry for file
            maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]]))
        except:
            pass
        if snpIndex % 100_000 == 0 and snpIndex != 0:
            print(snp, oldNucAtCurrentPos)
            print("time", time.time()-t1)
            print("index", snpIndex)
            t1 = time.time()
            print("progress",snpIndex / snpsLength)

    indexFile.write(str(snpPos) + "\n")


# cover for last case
for fileWithMissingSNP in filesWithoutSNP:
    dataToWrite[fileWithMissingSNP] += oldNucAtCurrentPos

print("printing dataToRight")

print(dataToWrite)
print("done with snps loop, now writing")
with open(outFilePath, "w") as outFile:
    for key in dataToWrite.keys():
        val = dataToWrite[key]
        outFile.write(">" + key.split("/")[-1] + "\n")
        outFile.write(val + "\n")
#
prof.disable()
prof.create_stats()
prof.print_stats()

# import pstats, io
# # from pstats import SortKey
# s = io.StringIO()
# # sortby = SortKey.CUMULATIVE
# ps = pstats.Stats(prof, stream=s)#.sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())
