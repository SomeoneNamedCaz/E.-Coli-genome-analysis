"""
NOTE: this isn't working yet
"""
from glob import glob
import time
import copy
import os
import sys
import cProfile
if len(sys.argv) < 3: # when you do args need ""
    print("requires a path to all of the vcf files like *.vcf and a output file name with path !use quotes!\n"
          "you may add a list of the types of snps to skip")
    exit(1)

if len(sys.argv) > 10:
    print("remember quotes")
    exit(1)

prof = cProfile.Profile()
prof.enable()

gsAlignPathString = sys.argv[1]#"./allGsAlignOutputs/*.vcf"
if gsAlignPathString.split(":")[0] == "MatchFile": # MatchFile:<file to match>:<pathPrefix>
    # this block is for using the same input files that were used in a previous run of this program
    # (if you say the file to match is the snp alignment output file)
    print("matching file")
    gsAlignFiles = []
    with open(gsAlignPathString.split(":")[1]) as file:
        for line in file:
            if line[0] == ">":
                gsAlignFiles.append(gsAlignPathString.split(":")[2] + line[1:].strip())
else:
    gsAlignFiles = glob(gsAlignPathString)


print(gsAlignFiles)


outFilePath = "/".join(sys.argv[2].split("/")[:-1]) + "/" #"./InsertAndDeleteCombinedGenomes/" # folder
outFileName = sys.argv[2].split("/")[-1]

try:
    thingsToSkip = sys.argv[3].split(",")
except IndexError:
    thingsToSkip = []


snps = [] # list of elements like this: [fileName, location, oldNuc, NewNuc, type]

numFiles = 0
for filePath in gsAlignFiles:
    numFiles += 1
    # if numFiles > 100:
    #     break
    print("opened one")
    with open(filePath) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or line[0] == "#": # or cols[7][5:] == "SUBSTITUTE":
                continue
            #               0                       1           2           3           4
            snps.append([filePath.split("/")[-1], int(cols[1]), cols[3], cols[4], cols[7][5:]])
            # if numFiles == 1:
            #     print([filePath, cols[1], cols[3], cols[4], cols[7][5:]])
            # snpPositions.add(cols[1])

# for snp in snps:
#     snp[0] = snp[0].split("/")[-1]

#     print(snp)
# print(len(snpPositions))

# with open("./SNPsList.tsv", "w") as file:
#     for snp in snps:
#         for val in snp:
#             file.write(val + "\t")
#         file.write("\n")


snps.sort(key=lambda a: a[1])
for snp in snps[:10000]:
    print(snp) #to test
# del snpPositions
allFiles = gsAlignFiles

for index in range(len(allFiles)):
    allFiles[index] = allFiles[index].split("/")[-1]
filesWithoutSNP = set(copy.deepcopy(allFiles))

removedFiles = set()
# print(filesWithoutSNP)
lastPos = snps[0][1]
# currFile = open(outFilePathPrefix + snps[0][1] + ".afa", "w")
oldNucAtCurrentPos = snps[0][2]
# longestNucLength = determineIfVariantIncludesNumOfNucleotideChange(currPos, 0) # -1 is a flag for no insertions or deletions
snpIndex = 0
dataToWrite = {} # {fileName:SNPlist}
print(len(snps))
snpsLength = len(snps) #* len(filesWithoutSNP)
for file in allFiles: # load the list with dicts with empty strings
    dataToWrite[file] = []#""
t1 = time.time()
numerrors = 0
numSkipped = 0

lastIndex = 0

# maxLenOfCurrSnp = 0 # max length of the snp that is currently in index to be added (i.e. first nucleotide has come and the last one hasn't)
maxLenOfSnpGenome = -1
snpIndex = -1
needToSkip = False
snpIndexPrintedToFile = False

if not os.path.exists(outFilePath):
    # print("tried to make dir")
    # exit(1)
    os.mkdir(outFilePath)

# add file name
outFilePath += outFileName # needs three letter extension
pathForSNPsIncludedIndexes = outFilePath[:-4] + "Indexes.txt"
lenBefore = 0
print(len(filesWithoutSNP))
print("start loop")
with open(pathForSNPsIncludedIndexes, "w") as indexFile, open(pathForSNPsIncludedIndexes[:-4]+"FrameShifted.txt", "w") as frameShiftIndexFile:
    for snp in snps: # list of elements like this: [fileName, location, oldNuc, NewNuc, type]
        snpIndex += 1
        if snp[4] in thingsToSkip:
            print(snp[4])
            continue
        snpPos = snp[1]
        if lastPos < snpPos: # if ran out of SNPs at the position
            # add the oldNuc for that SNP
            if not needToSkip: #419274
                # lenBefore = len(dataToWrite[filesWithoutSNP[0]]) # all should have same length
                for fileWithMissingSNP in filesWithoutSNP:
                    for char in oldNucAtCurrentPos:
                        dataToWrite[fileWithMissingSNP].append(char)

                # prob only need to do for the last one #TODO:check
                if len(filesWithoutSNP) != 0:
                    maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[fileWithMissingSNP]))

                for file in allFiles:
                    for char in "-" * (maxLenOfSnpGenome - len(dataToWrite[file])):
                        dataToWrite[file].append(char)
                for i in range(maxLenOfSnpGenome - lenBefore):
                    indexFile.write(str(lastPos + i/1000) + "\n")

                # # TODO: add this data. Could be first instead of last
                # if abs(len(oldNucAtCurrentPos) - len(dataToWrite[file])) % 3 != 0: # != 0
                #     frameShiftIndexFile.write(str(currPos) + "\n")
            try:
                if snpPos < snps[snpIndex + 9][1]: # if less than 10 snp at current position
                    needToSkip = True
                    continue
            except:
                pass
            # reset everything for next snp
            lastPos = snpPos
            filesWithoutSNP = filesWithoutSNP.union(removedFiles)  # add back removed files
            removedFiles = set()
            oldNucAtCurrentPos = snp[2]
            snpIndexPrintedToFile = False
            # get a random set element
            for randFile in filesWithoutSNP:
                break
            lenBefore = len(dataToWrite[randFile]) # should all have same length
        needToSkip = False
        try:
            filesWithoutSNP.remove(snp[0]) # throws error if trying to add duplicate snp
            removedFiles.add(snp[0])
            if len(snp[2]) > len(oldNucAtCurrentPos):
                oldNucAtCurrentPos = snp[2]
            maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]]) + len(snp[2]))
            if not snpIndexPrintedToFile and (len(snp[3]) - len(snp[2])) % 3 != 0:  # != 0
                frameShiftIndexFile.write(str(snpPos) + "\n")
                snpIndexPrintedToFile = True

            # add the snp nucs
            for char in snp[3]:
                dataToWrite[snp[0]].append(char) # add snp to entry for file
            # print("added this",snp[3])
            maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]]))
        except KeyError:
            pass
        if snpIndex % 2_000_000 == 0 and snpIndex != 0:
            print("time", time.time()-t1)
            print("index", snpIndex)
            t1 = time.time()
            print("progress",snpIndex / snpsLength)

    # cover for last case
    for fileWithMissingSNP in filesWithoutSNP:
        for char in oldNucAtCurrentPos:
            dataToWrite[fileWithMissingSNP].append(char)
        maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[fileWithMissingSNP]))
    for file in allFiles:
        for char in "-" * (maxLenOfSnpGenome - len(dataToWrite[file])):
            dataToWrite[file].append(char)
    for i in range(maxLenOfSnpGenome - lenBefore):
        indexFile.write(str(lastPos + i / 1000) + "\n")


print("printing dataToRight")

# print(dataToWrite)
print("done with snps loop, now writing")
with open(outFilePath, "w") as outFile:
    for key in dataToWrite.keys():
        val = "".join(dataToWrite[key])
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
