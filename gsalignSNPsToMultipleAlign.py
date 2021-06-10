"""
NOTE: this isn't working yet
"""
from glob import glob
import time
import copy
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
    # if numFiles == 4:
    #     break
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
filesWithoutSNP = glob(gsAlignPaths)
indexesWithoutSNPs = range(len(filesWithoutSNP))
for file in filesWithoutSNP:
    file = file.split("/")[-1]
filesWithoutSNPBackUp = copy.deepcopy(filesWithoutSNP)
# print(filesWithoutSNP)
currPos = snps[0][1]
outFilePathPrefix = "./substSNPCombinedGenomes/SNPMultipleAlignPos"
# currFile = open(outFilePathPrefix + snps[0][1] + ".afa", "w")
oldNucAtCurrentPos = snps[0][2]
# longestNucLength = determineIfVariantIncludesNumOfNucleotideChange(currPos, 0) # -1 is a flag for no insertions or deletions
snpIndex = 0
dataToWrite = {} # {fileName:SNPlist}
print(len(snps))
snpsLength = len(snps) #* len(filesWithoutSNP)
for file in filesWithoutSNP: # load the list with dicts with empty strings
    dataToWrite[file] = ""
t1 = time.time()
numerrors = 0
numSkipped = 0

lastIndex = 0
# for snp in snps:
#     snpPos = snp[1]
#     if currPos < snpPos:
#         print("snpIndex",snpIndex)
#         print("size of snps portions", snpIndex - lastIndex)
#         lastIndex = snpIndex
#         currPos = snpPos
#     snpIndex += 1
#
# exit(0)
snpIndex = -1
for snp in snps:# list of elements like this: [fileName, location, oldNuc, NewNuc, type]
    snpIndex += 1
    snpPos = snp[1]
    # if snpIndex % 100 == 0:
    #     print(snpIndex)
    if snp[4] != "SUBSTITUTE":
        # numSkipped += 1
        # if (numSkipped % 100 == 0):
        #     print("numskipped",numSkipped)
        continue
    if currPos < snpPos: # if ran out of SNPs at the position
        # print(snpPos-currPos)
        # add the oldNuc for that SNP
        for fileWithMissingSNP in filesWithoutSNP:
            dataToWrite[fileWithMissingSNP] += oldNucAtCurrentPos
            # snpIndex += 1
            # if snpIndex % 10_000_000 == 0:
            #     print(time.time() - t1)
            #     t1 = time.time()
            #     print("index", snpIndex)
            #     print(snpIndex / snpsLength)
        # reset everything for next snp
        currPos = snpPos
        filesWithoutSNP = copy.deepcopy(filesWithoutSNPBackUp)
        oldNucAtCurrentPos = snp[2]
        try:
            if snpPos < snps[snpIndex + 9][1]: # if less than 10 snp at current position
                continue
        except:
            0
    # fileName = snp[0]
    # if filesWithoutSNP.count(fileName) == 1:
    try:
        filesWithoutSNP.remove(snp[0]) # throws error if trying to add duplicate snp
        dataToWrite[snp[0]] += snp[3] # add snp to entry for file
    except:
        0
    #     numerrors += 1
    #     print("numerrors", numerrors)
    if snpIndex % 100_000 == 0 and snpIndex != 0:
        print("time", time.time()-t1)
        print("index", snpIndex)
        t1 = time.time()
        print("progress",snpIndex / snpsLength)
        # break
print("printing dataToRight")
print(dataToWrite)
print("done with snps loop, now writing")
with open("./substSNPcombinedGenomes/allSubstSNPsMoreThan9FastestMethod.afa", "w") as outFile:
    for key in dataToWrite.keys():
        val = dataToWrite[key]
        outFile.write(">" + key.split("/")[-1] + "\n")
        outFile.write(val + "\n")
#
prof.disable()
prof.create_stats()
print(prof.print_stats())

import pstats, io
# from pstats import SortKey
s = io.StringIO()
# sortby = SortKey.CUMULATIVE
ps = pstats.Stats(prof, stream=s)#.sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
