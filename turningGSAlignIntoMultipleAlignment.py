from functions import *
from copy import deepcopy

gsAlignPaths = "./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/gsAlignOutputs/*.maf"
firstFileAlignment = [] # [[refSection1, qry section1],[refsection2, qry section2], ...]
currentFileAlignment = []
# multipleGenomeAlignmentFilePath = "./multipleGenomAlignmentFromGSaligns.fasta"
multipleAlignment = [] # [[refSection1, qry1, qry2, qry3, ... ], ... ] in same order as first File segments
onFirstFile = True
# with open(multipleGenomeAlignmentFilePath, "w") as multipleAlignFile:
fileNumber = 1
numFilesToAlign = 3

# O(n^3) I think

def bruteForceSearch(currentFileAlignment, firstFileAlignment):
    possibleNewStartPointsInFirstFile = []
    possibleNewStartPointsInOtherFile = []
    numFound = 0
    for refSeg, qrySeg in currentFileAlignment:
        indexOfRefSegment = 0
        for refSeg2, qrySeg2 in firstFileAlignment:
            smallerLength = min(len(refSeg), len(refSeg2))
            # if not genesAreWithinPercentIdentical(refSeg[:20], refSeg2[:20]):  # this assumes that the segments are going to be shorter than 20bp
            #     continue

            if genesAreWithinPercentIdentical(refSeg[:20], refSeg2[:20]) and genesAreWithinPercentIdentical(
                    refSeg[:smallerLength], refSeg2[:smallerLength]):
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
                    smallerLength = min(smallerLength, len(multipleAlignment[indexOfRefSegment][0]))
                    # TODO: shift every other segment by deletions here
                    #  might be done
                    for index in range(smallerLength):
                        # checking for deletions
                        # print(multipleAlignment[indexOfRefSegment])
                        if refSeg[index] == "-" and multipleAlignment[indexOfRefSegment][0][  # first seq is reference
                            index] != "-":  # fails on case where aa segment in the multiple alignment is shorter than the refSeg
                            for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment])):
                                multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = insertAtIndex(
                                    multipleAlignment[indexOfRefSegment][indexOfAddedSegment], index, "-")
                        # shifting the qry seg over for - in multi align ref seg
                        if refSeg[index] != "-" and multipleAlignment[indexOfRefSegment][0][index] == "-":
                            qrySeg = insertAtIndex(qrySeg, index, "-")
                    for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment])):
                        # resize everything if it is smaller
                        multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = \
                            multipleAlignment[indexOfRefSegment][indexOfAddedSegment][:smallerLength]
                    multipleAlignment[indexOfRefSegment].append(qrySeg[:smallerLength])
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

# acutlally a newer bruteforce
def binarySearch(currentFileAlignment, firstFileAlignment):
    possibleNewStartPointsInFirstFile = []
    possibleNewStartPointsInOtherFile = []
    numFound = 0
    outerIndex = -1
    for refSeg_o, qrySeg_o in deepcopy(currentFileAlignment):
        outerIndex += 1
        refSeg = deepcopy(refSeg_o)
        qrySeg = deepcopy(qrySeg_o)
        indexOfRefSegment = 0

        for refSeg2_o, qrySeg2_o in deepcopy(firstFileAlignment): # copy so it aginmetm doesnt break things
            refSeg2 = deepcopy(refSeg2_o)
            qrySeg2 = deepcopy(qrySeg2_o)
            smallerLength = min(len(refSeg), len(refSeg2))
            if genesAreWithinPercentIdentical(refSeg[:20], refSeg2[:20]):  # this assumes that the segments are going to be shorter than 20bp
                #     continue
                for index in range(smallerLength - 1, -1, -1):  # could be a problem, note might have fixed it
                    if refSeg2[index] == "-" and refSeg[index] != "-":
                        refSeg = insertAtIndex(refSeg, index, "-")
                        qrySeg = insertAtIndex(qrySeg, index, "-")
                    if refSeg2[index] != "-" and refSeg[index] == "-":
                        refSeg2 = insertAtIndex(refSeg2, index, "-")
                        qrySeg2 = insertAtIndex(qrySeg2, index, "-")
                if  genesAreWithinPercentIdentical(refSeg[:smallerLength], refSeg2[:smallerLength]):
                    numFound += 1
                    if genesAreWithinPercentIdentical("TTTGACCTG", refSeg[:9]):
                        print("found")
                    # if (genesAreWithinPercentIdentical(refSeg[:smallerLength], "TTTGACCTGTCAGTAAAAGACACGTTGAAAACCTGGGGTCTGCTGGAGCTGGTCAACTCCGTGGTTGGGCTGATTATTGTGCTGATTATTAGCATGGTAGCGTAAAAGGACACCAGAGCCTGCCAATGGCAGGCTCAGACTGATGACAAAGTCTAAAAATGTGCCCGGACAGTCCCCTCGCCCCTCCGGGGAGAGGGTTAGGGTGAGGGGAACAGACCAGCAAAGTGCGCTTTGTTTATGCCAGATGCGACGCTGACGCGTCTTATCTGGCCTACAAAGGGCTAACGTGCAGGTTTTTAGCTTCAGGTAATATTGCGTACCAGCATTAGCAATGTCCTGTGATTTCTTTATTGATAAACAAAAGTCACGCCAATAATCGATTGCACATTCCCTGCAGTCACCTGCCCTCCGGTACGTGCATAATTTGCCGTTAATCCCAGACTCACCGCCGAAGTCCCTACTGCTCCTAACGATACCGTGTTATTCGCTGGAATAATCGTACCGTTGCGCGTCAACTGTACGCCGACGCCCTGTGCAGGTGAAAACGACGCGGTATTGGTGAAAATCGAGTTGCCCGCATCTGCGGTTGTGCCGGAGAGGTAATACCCCAGGTTTTGGCTTTTCGCACAATAAACGGTAAGAGGAATTGGCACTGAACCAGGGTAGTCCGGCAGAGTAACGGTGACATCACGAGCAGAAACATCGCAGCCGCCAGTAGGCACCACCACATCATTATTGGCGTAAATATTCCACACAAACTGGAAATCATCGCTGTTATAGTTGTTGGTCTGTCGCAAAATAAGCACGGCAATTAATGAGCCAGCTTTAATCGCCACCCCGCCCGCACTGCTCACAGGCGTCAAATAAAGCGCCACCGGCCACGGCTTATCCGTTCTCGAATTATAAACAACGCGCGGCGTTTCGCTGGTGGTAGGAAATGGATAGCTACTGCCACTATATTTTACGGTCCCGGAAAAATTAGATAACACGCCGCCATAAGCCGAGCCTCGTTGCAGTGTGACATAGTCTGTAATGGTTTCCGGATAATCGTTATGGCAAAAGATTTGCGTCGAAAGATCCACGACCAGGTTTTGCCCCACATTCACGACGGGCGCAAGGTTTACATAAACATTGGCGCTGCCACCGCCAATAGGGATAGCGGTACCATTGGCGGTTTTACAGGCGAATGACCAGGCATTTACCGACCAGCCCATCAGCAGTACAGCAAACAGGGTAATAACTCGTTTCATTACAATCATCTCTTCGGGTTCAGCTGTAGGTATAGGTGATGCTAATCACTGCCTGAATGGTTCCCTGAGTGGCTCCGCCATTTACTGTCAATGCTCTGACCTGTAACGGGAAGTGCGCTGATTGTGAGGAATCATCCACCTGAACTGTTTTGGTTGCGCCAGTATTCAATGTGTTGCCACTGTCATCCTGTAGCTCTAACTGGATGTTTTGCGCGGTCCCCTGGTTTTTATAATATCCGGTACTGTCGGCTGCCCCGCTGAAGCTGGCAGTGACCCTCGACGTTCCCACCGGACAATTAGTCAACTCAAGCGCAACATCATGCCAGGCCGATGCCGCCCCGGCAGACATAAGACTGAAAGAATAAAGATCGCCGAGATCAACCGTGGCATTGGTGGTGGAAACCGTACACGGTTTGGCGACGACCTTACCGTTCACCGTGATGGTGACATCGGCTGCCTGTATCGTCGCACTTGCGAGCGCCAATATTGCCGCCAATACATACCCACGTTTGCACCATTTCATGAGCATCTCCAGTTACTGATATTCAAGAGTGAAGGTAGCCGTGGCATTGATATGCCCCGCAGTGACAGGCACCTGTGTCGCCATTAGCCGGGCGTAAAAATTCAGCGTATTTGGTTTACCCGGCGTCAGGGTCGTCCACGAAAGCGCGGACGATGGAGCATTAAGGGGTATTTGATTTTGCTGCTCATTCAGAAGCTGTATTCCCAGTCCCGAAGCCGCTGACACCGTATTTTCAAGTGCAAGCAGGTTGGCATTGTGGCTATCTGCAACGCCAGTAAACCCAACCTTTACGGCAGAAACGGCATTACCACAGGGTGACAGCAAAATACGAAATGGAACAACAGGAGTCGTCGCGCCAATGTTGTTAAATTGCTTCGCCGCGTTTTCCATCAGATCAACAGTAAAATTGGTTGATTCAGCGGCCACACTACAGCCGTTATCCCTGACATAGCCGCGGATAGTAATCGTGCTATCCGCAGCCAAAGCGTGACTCACCGCCAGCCACAAAAAAGCGCACAGAAGATAAAAAGGTTTGTTTCTCATCACGCCCCCTTAACGACATTCAGCTGATAGCTGGGTTAATAACTGCTGCTGACTCTCTGGTGGCAGTTGATAATTGGCGACACAGTGAGCATTTTCCTCTTCTCCCCATTTCACCTGAACTTTTCCCGCTAAAGGCATTCCGCTGAGGTAAACCTGACCATTATCCGCAACAATGCCGCTACTCTGGCTACTCTCTGATGTCACCATCGCCCCAAACGGCAGCGGCTTATTATTGTGGGTCAGCGTCATGAGCAGTTTTATCCCAACGCGCGCTTTAAACTCTGCTCGCACGATCGCCCCACGAGTGGGAACAACGTTAGCAACCGCGTTATCTAAATCGACGTTATCAGCCAGGGTATTGGTATCCAGCGCCACTCTATTTTCCCGATATTCAGTGGCATAAGGCAGCACGGCATAACCACGCCAGTCGGTACGCACCCCCGTCTGGTTTTCGACTTTTGCATCTTTTGCGCCAGGCGCTTTAACAAGCACCACCGTATCGTTTAACGGCTGCCCCAGCGTTACGCCATTGGCATGAGCCAGTACCCCACCGCTGACTCCGTAATAGAGCTGCTTAATATCATCGCTATGGCTGTA"[:smallerLength])):
                    #     print("problem")
                    if len(multipleAlignment[indexOfRefSegment]) == 0:
                        multipleAlignment[indexOfRefSegment] = [refSeg[:smallerLength], qrySeg2[:smallerLength],
                                                                qrySeg[:smallerLength]]
                    else:
                        smallerLength = min(smallerLength, len(multipleAlignment[indexOfRefSegment][0]))
                        # TODO: shift every other segment by deletions here
                        #  might be done
                        for index in range(smallerLength):
                            # checking for deletions
                            # print(multipleAlignment[indexOfRefSegment])
                            if refSeg[index] == "-" and multipleAlignment[indexOfRefSegment][0][ # first seq is reference
                                index] != "-": # check for insertions
                                for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment]) - 1, -1, -1):
                                    multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = insertAtIndex(
                                        multipleAlignment[indexOfRefSegment][indexOfAddedSegment], index, "-")
                            # shifting the qry seg over for - in multi align ref seg
                            if refSeg[index] != "-" and multipleAlignment[indexOfRefSegment][0][index] == "-":
                                qrySeg = insertAtIndex(qrySeg, index, "-")
                        for indexOfAddedSegment in range(len(multipleAlignment[indexOfRefSegment])):
                            # resize everything if it is smaller
                            multipleAlignment[indexOfRefSegment][indexOfAddedSegment] = \
                            multipleAlignment[indexOfRefSegment][indexOfAddedSegment][:smallerLength]
                        multipleAlignment[indexOfRefSegment].append(qrySeg[:smallerLength]) # added twice on some segments added once on others (seems to happen when the bottom two (before the extra addition) segments match
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
    currentFileAlignment = []
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
                    multipleAlignment.append([]) # load multiple alignment with the right number for the first genom
                else:
                    currentFileAlignment.append([refSection, qrySection])
                    # print(refSection)
                    # find matching segments segments
                    # handle insertions
    ## after close of file
    firstFileAlignment.sort(key=lambda a: a[0])
    # print(len(firstFileAlignment))
    binarySearch(currentFileAlignment, firstFileAlignment)
    lensOfSegs = []
    for a in multipleAlignment:
        if a != []:
            lensOfSegs.append(len(a))
    if len(lensOfSegs) != 0:
        print("avg lens of segs", sum(lensOfSegs) / len(lensOfSegs))
        print("min", min(lensOfSegs))
        print("max", max(lensOfSegs))
    else:
        print("empty")

    onFirstFile = False

numSeg = 0
lensOfSegs = []
for a in multipleAlignment:
    print("_______")
    if a != []:
        numSeg += 1
        lensOfSegs.append(len(a))
    for s in a:
        # time.sleep(0.1)
        print(s)

print("num aligned segments",numSeg)
print("avg lens of segs", sum(lensOfSegs)/len(lensOfSegs))
print("min", min(lensOfSegs))
print("max", max(lensOfSegs))
