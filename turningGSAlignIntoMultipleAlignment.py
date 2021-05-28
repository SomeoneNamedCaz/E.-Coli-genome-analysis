from functions import *

gsAlignPaths = "./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.maf"
firstFileAlignment = [] # [[refSection1, qry section1],[refsection2, qry section2], ...]
currentFileAlignment = []
# multipleGenomeAlignmentFilePath = "./multipleGenomAlignmentFromGSaligns.fasta"
multipleAlignment = [] # [[refSection1, qry1, qry2, qry3, ... ], ... ] in same order as first File segments
onFirstFile = True
# with open(multipleGenomeAlignmentFilePath, "w") as multipleAlignFile:
fileNumber = 1
numFilesToAlign = 3000

def bruteForceSearch(currentFileAlignment, firstFileAlignment):
    possibleNewStartPointsInFirstFile = []
    possibleNewStartPointsInOtherFile = []
    numFound = 0
    for refSeg, qrySeg in currentFileAlignment:
        indexOfRefSegment = 0
        for refSeg2, qrySeg2 in firstFileAlignment:
            smallerLength = min(len(refSeg), len(refSeg2))
            if not genesAreWithinPercentIdentical(refSeg[:20], refSeg2[
                                                               :20]):  # this assumes that the segments are going to be shorter than 20bp
                continue
            if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                numFound += 1
                if len(multipleAlignment[indexOfRefSegment]) == 0:
                    for index in range(smallerLength):
                        if refSeg2[index] == "-" and refSeg[index] != "-":
                            refSeg = insertAtIndex(refSeg, index, "-")
                            qrySeg = insertAtIndex(qrySeg, index, "-")
                        if refSeg2[index] != "-" and refSeg[index] == "-":
                            refSeg2 = insertAtIndex(refSeg2, index, "-")
                            qrySeg2 = insertAtIndex(qrySeg2, index, "-")
                    multipleAlignment[indexOfRefSegment] = [refSeg[:smallerLength], qrySeg2[:smallerLength],
                                                            qrySeg[:smallerLength]]
                else:
                    # TODO: fix out of range error

                    #                        if new seg or ref seg is smaller
                    for index in range(min(smallerLength, multipleAlignment[indexOfRefSegment][0])):
                        if refSeg[index] == "-" and multipleAlignment[indexOfRefSegment][0][  # first seq is reference
                            index] != "-":
                            for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment])):
                                multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = insertAtIndex(
                                    multipleAlignment[indexOfRefSegment][indexOfAddedSegment], index, "-")

                    multipleAlignment[indexOfRefSegment] += qrySeg
                if len(refSeg) > len(refSeg2):
                    possibleNewStartPointsInFirstFile.append(refSeg[smallerLength:])
                elif len(refSeg) < len(refSeg2):
                    possibleNewStartPointsInOtherFile.append(refSeg2[smallerLength:])
            indexOfRefSegment += 1

    print("second loop")  # these two loops didn't find anything
    for refSeg in currentFileAlignment:
        for refSeg2 in possibleNewStartPointsInFirstFile:
            smallerLength = min(len(refSeg), len(refSeg2))
            if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                print("found")
                numFound += 1
    for refSeg in possibleNewStartPointsInOtherFile:
        for refSeg2 in firstFileAlignment:
            smallerLength = min(len(refSeg), len(refSeg2))
            if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                # print("found")
                numFound += 1
    if currentFileAlignment != []:
        print(numFound)
        print(numFound / len(currentFileAlignment))

def binarySearch(currentFileAlignment, firstFileAlignment):
    possibleNewStartPointsInFirstFile = []
    possibleNewStartPointsInOtherFile = []
    numFound = 0
    for refSeg, qrySeg in currentFileAlignment:
        indexOfRefSegment = 0

        for refSeg2, qrySeg2 in firstFileAlignment:
            smallerLength = min(len(refSeg), len(refSeg2))
            # if not genesAreWithinPercentIdentical(refSeg[:20], refSeg2[:20]):  # this assumes that the segments are going to be shorter than 20bp
            #     continue
            if genesAreWithinPercentIdentical(refSeg[:20], refSeg2[:20]) and genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                numFound += 1
                if len(multipleAlignment[indexOfRefSegment]) == 0:
                    for index in range(smallerLength):
                        if refSeg2[index] == "-" and refSeg[index] != "-":
                            refSeg = insertAtIndex(refSeg, index, "-")
                            qrySeg = insertAtIndex(qrySeg, index, "-")
                        if refSeg2[index] != "-" and refSeg[index] == "-":
                            refSeg2 = insertAtIndex(refSeg2, index, "-")
                            qrySeg2 = insertAtIndex(qrySeg2, index, "-")
                    multipleAlignment[indexOfRefSegment] = [refSeg[:smallerLength], qrySeg2[:smallerLength],
                                                            qrySeg[:smallerLength]]
                else:
                    # TODO: shift every other segment by deletions here
                    #  might be done
                    for index in range(smallerLength):
                        if refSeg[index] == "-" and multipleAlignment[indexOfRefSegment][0][ # first seq is reference
                            index] != "-":
                            for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment])):
                                multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = insertAtIndex(
                                    multipleAlignment[indexOfRefSegment][indexOfAddedSegment], index, "-")

                    multipleAlignment[indexOfRefSegment].append(qrySeg)
                if len(refSeg) > len(refSeg2):
                    possibleNewStartPointsInFirstFile.append(refSeg[smallerLength:])
                elif len(refSeg) < len(refSeg2):
                    possibleNewStartPointsInOtherFile.append(refSeg2[smallerLength:])
            indexOfRefSegment += 1

    # print("second loop")  # these two loops didn't find anything
    # for refSeg in currentFileAlignment:
    #     for refSeg2 in possibleNewStartPointsInFirstFile:
    #         smallerLength = min(len(refSeg), len(refSeg2))
    #         if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
    #             print("found")
    #             numFound += 1
    # for refSeg in possibleNewStartPointsInOtherFile:
    #     for refSeg2 in firstFileAlignment:
    #         smallerLength = min(len(refSeg), len(refSeg2))
    #         if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
    #             # print("found")
    #             numFound += 1
    # if currentFileAlignment != []:
    #     print(numFound)
    #     print(numFound / len(currentFileAlignment))


for path in glob(gsAlignPaths):
    if fileNumber > numFilesToAlign:
        break
    print("on file", fileNumber, "file path", path)
    fileNumber += 1
    # load files
    with open(path) as file:
        refSection = ""
        qrySection = ""
        for line in file:
            line = line.strip()
            if line == "":
                continue

            if line[:5] == "s ref":
                # in reference seq
                refSection = line.split()[-1]
            elif line[:5] == "s qry":
                # in query seq
                qrySection = line.split()[-1]
                if onFirstFile:
                    firstFileAlignment.append([refSection, qrySection])
                    multipleAlignment.append([]) # load multiple alignment with the right number for the first genome
                else:
                    currentFileAlignment.append([refSection, qrySection])
                    # print(refSection)
                    # find matching segments segments
                    # handle insertions
    ## after close of file
    firstFileAlignment.sort(key=lambda a: a[0])
    # print(len(firstFileAlignment))
    binarySearch(currentFileAlignment, firstFileAlignment)


    onFirstFile = False

for a in multipleAlignment:
    print("_______")
    for s in a:
        time.sleep(0.1)
        print(s)

