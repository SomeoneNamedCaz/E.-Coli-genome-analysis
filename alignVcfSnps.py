
from glob import glob
import time
import copy
import os
import sys
import cProfile
from secondaryPythonScripts.functions import *


def fixNucsOnOtherStrand(refGenomeSeq, location, refNucs, altNucs, snpType):
    correspondingRefSeq = refGenomeSeq[location:len(refNucs) + location]
    if refGenomeSeq[location:len(refNucs) + location] != refNucs:
        # if snpType == "INSERT":
        #     location -= len(altNucs) - 1
        if snpType == "DELETE":
            location -= len(refNucs) - 1
        if refNucs.count("N") == 0 and altNucs.count("N") == 0:
            refNucs = reverseComplement(refNucs)
            altNucs = reverseComplement(altNucs)
            if refGenomeSeq[location:len(refNucs) + location] != refNucs:
                print("real problem")
        else:
            print("problem", refNucs, altNucs)
    return location, refNucs, altNucs

def readInSnps(listOfFiles, refGenomeSeq, ignoreRefSeq):
    numFiles = 0
    snps = []  # list of elements like this: [fileName, location, oldNuc, NewNuc, type]
    for filePath in listOfFiles:
        numFiles += 1
        # if numFiles > 100:
        #     break
        with open(filePath) as file:
            for line in file:
                line = line.strip()
                cols = line.split("\t")
                if line == "" or line[0] == "#":  # or cols[7][5:] == "SUBSTITUTE":
                    continue
                
                location = int(cols[1]) - 1
                refNucs = cols[3]
                altNucs = cols[4]
                if refNucs.count("N") != 0 or altNucs.count("N") != 0:
                    continue
                snpType = cols[7][5:]
                if not ignoreRefSeq:
                    location, refNucs, altNucs = fixNucsOnOtherStrand(refGenomeSeq, location, refNucs, altNucs, snpType)
                
                #               0                       1            2      3         4
                snps.append([filePath.split("/")[-1], location, refNucs, altNucs, snpType])
    return snps
def alignVcfSnps(gsAlignPathString, outFilePath, thingsToSkip=[], numSnpsRequired=10, debug = False, refSeqPath="/Users/cazcullimore/Documents/ericksonLabCode/refGenomes/k-12.fasta", ignoreRefSeq=True, includeGenomeWithoutAnySnps=True):

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

    if debug:
        print(gsAlignFiles)

    outFileName = outFilePath.split("/")[-1]

    refSeq = readInFastaAsList(refSeqPath)[1]

    outFileDirectory = "/".join(outFilePath.split("/")[:-1]) + "/"
    ## parse files
    snps = readInSnps(gsAlignFiles, refSeq, ignoreRefSeq=ignoreRefSeq)
    snps.sort(key=lambda a: a[1])
    # if debug:
    for snp in snps[:1_00]:
        print(snp) #to test

    allFiles = gsAlignFiles

    for index in range(len(allFiles)):
        allFiles[index] = allFiles[index].split("/")[-1]
    allFiles = set(allFiles)
    if includeGenomeWithoutAnySnps:
        allFiles.add("genomeWitoutAnySnps")
    filesWithoutSNP = set(copy.deepcopy(allFiles))

    removedFiles = set()
    # print(filesWithoutSNP)
    lastPos = snps[0][1]
    # currFile = open(outFilePathPrefix + snps[0][1] + ".afa", "w")
    oldNucsAtCurrentPos = snps[0][2]
    # longestNucLength = determineIfVariantIncludesNumOfNucleotideChange(currPos, 0) # -1 is a flag for no insertions or deletions
    snpIndex = 0
    dataToWrite = {} # {fileName:SNPlist}
    if debug:
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

    if not os.path.exists(outFileDirectory):
        # print("tried to make dir")
        # exit(1)
        os.mkdir(outFileDirectory)


    def alignSnps(lengthOfLongestSnpGenome, longestRefNucSeq, filesToSnpGenome, filesWithoutSNP, lenBefore, lastPosition, indexFile, fileNameToLengthOfInsert, lengthOfLongestInsert,numberOfNucelotidesToStopAlignmentAt=-1):
        if debug:
            testPositions = [1708844, 1708847]
            if lastPosition in testPositions:
                print("error")
                for k, v in dataToWrite.items():
                    print(k[5:15] + "\t" + str(len(v)) + "\t"+ "".join(v[-100:]))
                print(longestRefNucSeq)
                
        for fileWithMissingSNP in filesWithoutSNP:
            if not len(filesToSnpGenome[fileWithMissingSNP]) > lenBefore: # in other words if equal to len befor
                filesToSnpGenome[fileWithMissingSNP].append(longestRefNucSeq[0])
        maxLenAfterDashes = lengthOfLongestSnpGenome

        
        if debug:
            if lastPosition in testPositions:
                print("error")
                for k, v in dataToWrite.items():
                    print(k[5:15] + "\t" + str(len(v)) + "\t"+ "".join(v[-100:]))

        numIndexesAdded = 0
        indexesAdded = []
        # add insert dashes
        # After this, all genomes will be the same length if there is no delete.
        # If there is a delete, it will be the longest genome that's why we can add reference nucleotides to balance the lengths
        for file in allFiles:
            # dictionary has the keys of all of the inserts
            if file in fileNameToLengthOfInsert.keys():
                # add dashes to potentially smaller inserts
                for _ in range(lengthOfLongestInsert - fileNameToLengthOfInsert[file]):
                    filesToSnpGenome[file].append("-")
            else:
                for _ in range(lengthOfLongestInsert):
                    filesToSnpGenome[file].append("-")
            if len(filesToSnpGenome[file]) > maxLenAfterDashes:
                maxLenAfterDashes = len(filesToSnpGenome[file])



        for i in range(lengthOfLongestInsert + 1): # add insert indexes + subst index (the plus 1)
            numIndexesAdded += 1
            indexFile.write(str(lastPosition + i / 1000) + "\n")
            indexesAdded.append(str(lastPosition + i / 1000) + "\n")

        lengthOfLongestSnpGenome = maxLenAfterDashes

        maxLen = numberOfNucelotidesToStopAlignmentAt + lenBefore + lengthOfLongestInsert
        numRefNucsAdded = 0
        for file in allFiles:
            if numberOfNucelotidesToStopAlignmentAt == -1:
                numRefNucsToAdd = lengthOfLongestSnpGenome - len(filesToSnpGenome[file])
            else:
                numRefNucsToAdd = min(lengthOfLongestSnpGenome, maxLen) - len(filesToSnpGenome[file])
            if numRefNucsToAdd <= 0: # if already too long or at proper length don't extend (fixes problems with sign changes)
                continue

            if numberOfNucelotidesToStopAlignmentAt == -1:
                for char in longestRefNucSeq[-numRefNucsToAdd:]:
                    filesToSnpGenome[file].append(char)
                if numRefNucsToAdd > numRefNucsAdded:
                    numRefNucsAdded = numRefNucsToAdd
            else:
                for char in longestRefNucSeq[:numberOfNucelotidesToStopAlignmentAt][-numRefNucsToAdd:]:
                    filesToSnpGenome[file].append(char)
                if numRefNucsToAdd > numRefNucsAdded:
                    numRefNucsAdded = numRefNucsToAdd

        for i in range(min(numRefNucsAdded,len(longestRefNucSeq))): # add delete indexes if there are any
            numIndexesAdded += 1
            indexFile.write(str(lastPosition + 1 + i) + "\n") # because of adding first ref nuc

            indexesAdded.append(str(lastPosition + 1 + i) + "\n")

        if debug:
            print(numRefNucsAdded, "numRefNucsAdded")
            print("longestRefNucSeq",longestRefNucSeq)


        if debug:
            if len(dataToWrite[file]) - lenBefore != numIndexesAdded and numberOfNucelotidesToStopAlignmentAt == -1:
                print(len(dataToWrite[file]) - lenBefore,"and",numIndexesAdded, "differ at", lastPosition)#, dataToWrite)
            if lastPosition in testPositions:
                print("numRefNucsAdded", numRefNucsAdded)
                print("after align")
                for k, v in dataToWrite.items():
                    print(k[5:15] + "\t" + str(len(v)) + "\t" + "".join(v[-100:]))# str(len(v)) should all be the same here
                print(indexesAdded, "done,", longestRefNucSeq[numRefNucsAdded+1:], "really done now")

        return longestRefNucSeq[numRefNucsAdded+1:] # plus one because of added the first nucleotide earlier




    # add file name
    outFileDirectory += outFileName # needs three letter extension
    pathForSNPsIncludedIndexes = outFileDirectory[:-4] + "Indexes.txt"
    lenBefore = 0
    overlappingOldNucs = ""
    fileNameToLengthOfInsert = {} # string : int
    lengthOfLongestInsert = 0
    fileNameToNextValidPosition = {} # only for conflicts within GSalign files
    for file in allFiles:
        fileNameToNextValidPosition[file] = 0
    if debug:
        print(len(filesWithoutSNP))
        print("start loop")


    ######################## MAIN LOOP ###########################
    with open(pathForSNPsIncludedIndexes, "w") as indexFile, open(pathForSNPsIncludedIndexes[:-4] + "FrameShifted.txt", "w") as frameShiftIndexFile:
        for snp in snps: # list of elements like this: [fileName, location, oldNuc, NewNuc, type]
            if debug:
                print(snp)
            snpIndex += 1
            snpType = snp[4]
            if snpType in thingsToSkip:
                print("skipped:", snp[4])
                continue
            snpPos = snp[1]
            refNucs = snp[2]
            altNucs = snp[3]
            # check for using the other strand
            if not ignoreRefSeq:
                if refSeq[snpPos:len(refNucs) + snpPos] != refNucs:
                    if len(refNucs) > 1 or len(altNucs) > 1:
                        if debug:
                            print("nucs on wrong strand longer than 1 nuc. snp data: " + str(snp))
                        continue # we need to move the position back, so we probably should do that earlier in the program
                                 # before we sort the snps
                        # raise Exception("nucs on wrong strand longer than 1 nuc. snp data: " + str(snp))
                    refNucs = reverseComplement(refNucs)
                    altNucs = reverseComplement(altNucs)
                    if refSeq[snpPos:len(refNucs) + snpPos] != refNucs:
                        print("broken")
                        continue

            if lastPos < snpPos: # if ran out of SNPs at the position
                # add the oldNuc for that SNP
                try:
                    if snpPos < snps[snpIndex + numSnpsRequired - 1][1]: # if less than 10 snp at current position
                        needToSkip = True
                        continue
                except:
                    pass

                numNucsToAlign = -1
                if lastPos + len(oldNucsAtCurrentPos) > snpPos: # if overlap is starting or continuing
                    numNucsToAlign = snpPos - lastPos
                oldNucsAtCurrentPos = alignSnps(maxLenOfSnpGenome, oldNucsAtCurrentPos, dataToWrite, filesWithoutSNP,
                                               lenBefore, lastPos, indexFile, fileNameToLengthOfInsert,
                                               lengthOfLongestInsert,numNucsToAlign)

                # reset everything for next snp position
                lastPos = snpPos
                filesWithoutSNP = filesWithoutSNP.union(removedFiles)  # add back removed files
                removedFiles = set()
                fileNameToLengthOfInsert = {}
                if len(snp[2]) > len(oldNucsAtCurrentPos):
                    oldNucsAtCurrentPos = snp[2]
                snpIndexPrintedToFile = False
                lengthOfLongestInsert = 0
                # get a random set element

                for randFile in filesWithoutSNP:
                    break
                lenBefore = len(dataToWrite[randFile]) # get smallest
                for file in filesWithoutSNP:
                    if lenBefore > len(dataToWrite[file]):
                        lenBefore = len(dataToWrite[file])
            needToSkip = False
            try:
                if fileNameToNextValidPosition[snp[0]] > snpPos:
                    if debug:
                        print("internal vcf conflict")
                    continue
                filesWithoutSNP.remove(snp[0]) # throws error if trying to add duplicate snp
                removedFiles.add(snp[0])
                if len(snp[2]) > len(oldNucsAtCurrentPos):
                    oldNucsAtCurrentPos = snp[2]
                maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]])) #+ len(snp[2]))
                if not snpIndexPrintedToFile and (len(snp[3]) - len(snp[2])) % 3 != 0:  # != 0
                    frameShiftIndexFile.write(str(snpPos) + "\n")
                    snpIndexPrintedToFile = True

                # add the snp nucs
                for char in snp[3]:
                    dataToWrite[snp[0]].append(char) # add snp to entry for file
                maxLenOfSnpGenome = max(maxLenOfSnpGenome, len(dataToWrite[snp[0]]))
                if snpType == "DELETE":
                    for i in range(len(refNucs) - len(altNucs)):
                        dataToWrite[snp[0]].append("-")
                elif snpType == "INSERT":
                    insertLength = len(altNucs) - len(refNucs)
                    fileNameToLengthOfInsert[snp[0]] = insertLength
                    if insertLength > lengthOfLongestInsert:
                        lengthOfLongestInsert = insertLength
                # print("added this",snp[3])
                fileNameToNextValidPosition[snp[0]] = len(refNucs) + snpPos
            except KeyError: # prob not needed anymore
                pass
            if snpIndex % 2_000_000 == 0 and snpIndex != 0:
                print("time", time.time()-t1)
                print("index", snpIndex)
                t1 = time.time()
                print("progress",snpIndex / snpsLength)

        # cover for last case
        # if lastPos + len(oldNucsAtCurrentPos) + len(overlappingOldNucs) > snpPos:  # if next snps are within the position
        #     numNucsToAlign = snpPos - lastPos
        # if overlappingOldNucs != "":
        #     overlappingOldNucs = alignSnps(maxLenOfSnpGenome, overlappingOldNucs + oldNucsAtCurrentPos, dataToWrite, filesWithoutSNP, lenBefore, lastPos, indexFile, numNucsToAlign)
        # alignSnps(maxLenOfSnpGenome, overlappingOldNucs, dataToWrite, filesWithoutSNP, lenBefore, lastPos + numNucsToAlign, indexFile)
        if debug:
            print(overlappingOldNucs)
        alignSnps(maxLenOfSnpGenome, oldNucsAtCurrentPos, dataToWrite, filesWithoutSNP, lenBefore, lastPos,
                  indexFile, fileNameToLengthOfInsert, lengthOfLongestInsert)
    # print("printing dataToRight")
    if debug:
        print(dataToWrite)
        print("done with snps loop, now writing")
    # if __name__ == "__main__":
    with open(outFileDirectory, "w") as outFile:
        for key in dataToWrite.keys():
            val = "".join(dataToWrite[key])
            outFile.write(">" + key.split("/")[-1] + "\n")
            outFile.write(val + "\n")
    return dataToWrite
if __name__ == "__main__":
    if len(sys.argv) < 3:  # when you do args need ""
        print("requires a path to all of the vcf files like *.vcf and a output file name with path !use quotes!\n"
              "you may add a list of the types of snps to skip and the cutoff to remove infrequent snps")
        exit(1)

    if len(sys.argv) > 10:
        print("remember quotes")
        exit(1)

    try:
        thingsToSkip = sys.argv[3].split(",")
    except IndexError:
        thingsToSkip = []

    try:
        numSnpsRequired = sys.argv[4]
    except IndexError:
        numSnpsRequired = 1

    gsAlignPath = sys.argv[1]  # "./allGsAlignOutputs/*.vcf"
    outFilePath = sys.argv[2]


    alignVcfSnps(gsAlignPath, outFilePath, thingsToSkip,numSnpsRequired,debug=False)

