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
    print("requires a path to all of the vcf files like *.vcf and a output file name/path")
    exit(1)

prof = cProfile.Profile()
prof.enable()

gsAlignPaths = sys.argv[1]#"./allGsAlignOutputs/*.vcf"

outFilePath = "/".join(sys.argv[2].split("/")[:-1]) + "/" #"./InsertAndDeleteCombinedGenomes/" # folder
outFileName = sys.argv[2].split("/")[-1]

snps = [] # list of elements like this: [fileName, location, oldNuc, NewNuc, type]

numFiles = 0
for filePath in glob(gsAlignPaths):
    numFiles += 1
    # if numFiles > 200:
    #     break
    # print("opened one")
    with open(filePath) as file:
        for line in file:
            line = line.strip()
            cols = line.split("\t")
            if line == "" or line[0] == "#": # or cols[7][5:] == "SUBSTITUTE":
                continue
            #               0                       1           2           3           4
            snps.append([filePath, int(cols[1]), cols[3], cols[4], cols[7][5:]])
            # if numFiles == 1:
            #     print([filePath, cols[1], cols[3], cols[4], cols[7][5:]])
            # snpPositions.add(cols[1])

for snp in snps:
    snp[0] = snp[0].split("/")[-1]

#     print(snp)
# print(len(snpPositions))

# with open("./SNPsList.tsv", "w") as file:
#     for snp in snps:
#         for val in snp:
#             file.write(val + "\t")
#         file.write("\n")

#FIXME:I don't thing this works since the added positions doesn't have a real position in the genome since they are inserts or deletes
# snpsWithInserstAsSingleNucs = []
# for snp in snps:
#     if snp[4] == "INSERT":
#         for nuc in snp[2]:
#             snpsWithInserstAsSingleNucs.append([])
#     if snp[4] == "DELETE":

snps.sort(key=lambda a: a[1])
# del snpPositions
allFiles = glob(gsAlignPaths)

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
print(snps[:1000])
snpsLength = len(snps) #* len(filesWithoutSNP)
for file in allFiles: # load the list with dicts with empty strings
    dataToWrite[file] = ""
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
    os.mkdir(outFilePath)

# add file name
outFilePath += outFileName # needs three letter extension
pathForSNPsIncludedIndexes = outFilePath[:-4] + "Indexes.txt"
lenBefore = 0
print("start loop")
with open(pathForSNPsIncludedIndexes, "w") as indexFile, open(pathForSNPsIncludedIndexes[:-4]+"FrameShifted.txt", "w") as frameShiftIndexFile:
    for snp in snps: # list of elements like this: [fileName, location, oldNuc, NewNuc, type]
        snpIndex += 1
        snpPos = snp[1]
        if lastPos < snpPos: # if ran out of SNPs at the position
            # add the oldNuc for that SNP
            if not needToSkip: #419274
                # lenBefore = len(dataToWrite[filesWithoutSNP[0]]) # all should have same length
                for fileWithMissingSNP in filesWithoutSNP:
                    dataToWrite[fileWithMissingSNP] += oldNucAtCurrentPos

                # prob only need to do for the last one #TODO:check
                maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[fileWithMissingSNP]))

                for file in allFiles:
                    dataToWrite[file] += "-" * (maxLenOfSnpGenome - len(dataToWrite[file]))
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
            filesWithoutSNP.union(removedFiles)  # add back removed files
            removedFiles = set()
            oldNucAtCurrentPos = snp[2]
            snpIndexPrintedToFile = False
            # get a random set element
            for randFile in filesWithoutSNP:
                break
            lenBefore = len(dataToWrite[randFile])
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
            dataToWrite[snp[0]] += snp[3] # add snp to entry for file
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
    indexFile.write(str(snpPos) + "\n")
    for fileWithMissingSNP in filesWithoutSNP:
        dataToWrite[fileWithMissingSNP] += oldNucAtCurrentPos
        maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[fileWithMissingSNP]))
    for file in allFiles:
        dataToWrite[file] += "-" * (maxLenOfSnpGenome - len(dataToWrite[file]))


print("printing dataToRight")

# print(dataToWrite)
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
