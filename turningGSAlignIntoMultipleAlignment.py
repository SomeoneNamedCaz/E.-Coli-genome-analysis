from functions import *

gsAlignPaths = "./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.maf"
firstFileAlignment = [] # [[refSection1, qry section1],[refsection2, qry section2], ...]
currentFileAlignment = []
# multipleGenomeAlignmentFilePath = "./multipleGenomAlignmentFromGSaligns.fasta"
multipleAlignment = [] # [[refSection1, qry1, qry2, qry3, ... ], ... ] in same order as first File segments
onFirstFile = True
# with open(multipleGenomeAlignmentFilePath, "w") as multipleAlignFile:
for path in glob(gsAlignPaths):

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
                    firstFileAlignment.append([re.sub("-","",refSection), qrySection])
                else:
                    currentFileAlignment.append([re.sub("-","",refSection), qrySection])
                    # find matching segments segments
                    # handle insertions


    ## after close of file

    numFound = 0

    possibleNewStartPointsInFirstFile = []
    possibleNewStartPointsInOtherFile = []
    sum = 0
    for refSeg, qrySeg in currentFileAlignment:
        sum += len(refSeg)
        indexOfRefSegment = 0
        for refSeg2, qrySeg2 in firstFileAlignment:
            smallerLength = min(len(refSeg), len(refSeg2))
            if genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                print("found")
                numFound += 1

                for index in refSeg2[:smallerLength]:
                    if refSeg2[index] == "-" and refSeg[index] != "-":
                        refSeg = insertAtIndex(refSeg, index, "-")
                        qrySeg = insertAtIndex(qrySeg, index, "-")
                    if refSeg2[index] != "-" and refSeg[index] == "-":
                        refSeg2 = insertAtIndex(refSeg2, index, "-")
                        qrySeg2 = insertAtIndex(qrySeg2, index, "-")
                if len(multipleAlignment[indexOfRefSegment]) == 0:
                    multipleAlignment.append([refSeg,qrySeg2,qrySeg])
                else:
                    multipleAlignment[indexOfRefSegment] += qrySeg
                    # TODO: shift every other segment by deletions here
                if len(refSeg) > len(refSeg2):
                    possibleNewStartPointsInFirstFile.append(refSeg[smallerLength:])
                elif len(refSeg) < len(refSeg2):
                    possibleNewStartPointsInOtherFile.append(refSeg2[smallerLength:])
            indexOfRefSegment += 1

    print(sum)
    print("second loop") # these two loops didn't find anything
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
                print("found")
                numFound += 1
    if currentFileAlignment != []:
        print(numFound)
        print(numFound / len(currentFileAlignment))
        break

    onFirstFile = False