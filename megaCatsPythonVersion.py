"""Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups	Results_Filename
 313 	 7.46929814283269 	 0.0238815503170304 	 2 	 Y 	 group1(19_A,_18_C,_1_T)|group2(8_A,_30_C,_1_T) 	severity-rMsaInput.txt-rResultChisqTest.txt"""
import sys
from scipy import stats
import numpy as np

alignedFilePath = sys.argv[1]
metadataFilePath = sys.argv[2]

class NucCounter():

    def __init__(self,nameOfGroup = "unnamed", counts=[]):
        self._nucCounts = counts
        self.nameOfGroup = nameOfGroup
    def getDistributionAt(self, positionInGenome):
        return self._nucCounts[positionInGenome]
    def appendDistribution(self, distribution):
        self._nucCounts.append(distribution)

possibleClassValues = [] # maybe do a separate class for this
seqNameToMetaDataType = {}
with open(metadataFilePath) as metadataFile: # process metadata
    isFirstLine = True
    for line in metadataFile:
        line = line.strip()
        cols = line.split("\t")
        if isFirstLine:
            for _ in range(len(cols)-1):
                possibleClassValues.append(set())
            isFirstLine = False
        else:
            seqNameToMetaDataType[cols[0]] = cols[1:]
            for colIndex in range(len(cols[1:])):
                # try:
                possibleClassValues[colIndex].add(cols[colIndex + 1])
                # except TypeError:
                #     possibleClassValues[colIndex] = {}

print(possibleClassValues)
nucCounters = {}
# maybe better with a matrix but 2 minutes isn't that bad

with open(alignedFilePath) as alignedFile:
    # for metaDataTypes in range(len(possibleClassValues)):
    #     nucCounters =
    # hard code for now

    for val1 in possibleClassValues[0]:
        nucCounters[val1] = NucCounter(nameOfGroup=val1)
        # for val2 in possibleClassValues[1]:
        #     nucCounters[val1+val2] = NucCounter(nameOfGroup=val1+val2)

            # # TODO: can have conflict if possibleClassValues[0] and possibleClassValues[1] have same values
    for val2 in possibleClassValues[1]:
        nucCounters[val2] = NucCounter(nameOfGroup=val2)
    nameOfSeq = "unknown"
    for line in alignedFile:
        line = line.strip()
        if line[0] == ">":
            nameOfSeq = line[1:]
        else:
            # nameOfNucCounter = "".join(seqNameToMetaDataType[nameOfSeq])

            for i in range(len(possibleClassValues)):
                nameOfNucCounter = seqNameToMetaDataType[nameOfSeq][i]
                charIndex = -1
                for char in line:
                    charIndex += 1
                    try:
                        nucCounters[nameOfNucCounter]._nucCounts[charIndex][char] += 1
                    except KeyError:
                        nucCounters[nameOfNucCounter]._nucCounts[charIndex][char] = 1
                    except IndexError:
                        nucCounters[nameOfNucCounter]._nucCounts.append({char:1})


for metaDataCategoryIndex in range(len(possibleClassValues)):
    nucCountersOfOneMetadataCatagory = []
    for option in possibleClassValues[metaDataCategoryIndex]:
        # ["A","T","C","G"]
        # for
        nucCountersOfOneMetadataCatagory.append(nucCounters[option])
    # TODO: only works with two options in each metadatacatagory
    nucIndex = 0
    while True:
        try:
            dist1 = nucCountersOfOneMetadataCatagory[0].getDistributionAt([nucIndex])
            dist2 = nucCountersOfOneMetadataCatagory[1].getDistributionAt([nucIndex])
            for value in range(max(len(dist1.values), len(dist2.values))):

        except IndexError:
            break