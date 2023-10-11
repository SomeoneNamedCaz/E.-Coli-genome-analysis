"""Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups	Results_Filename
 313 	 7.46929814283269 	 0.0238815503170304 	 2 	 Y 	 group1(19_A,_18_C,_1_T)|group2(8_A,_30_C,_1_T) 	severity-rMsaInput.txt-rResultChisqTest.txt"""
import sys
from scipy import stats
import numpy as np

alignedFilePath = sys.argv[1]
metadataFilePath = sys.argv[2]

charToIndex = {"A":0, "T":1, "C":2, "G":3, "-": 4, "N":5}
class NucCounter():

    def __init__(self,countLength = 0, nameOfGroup = "unnamed",):
        self.nucCounts = []
        for _ in range(countLength):
            self.nucCounts.append([0] * len(charToIndex.values()))
        self.nameOfGroup = nameOfGroup
    def getDistributionAt(self, positionInGenome):
        return self.nucCounts[positionInGenome]
    def appendDistribution(self, distribution):
        self.nucCounts.append(distribution)

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

def makeNucCounters(lengthOfSeqs, optionsForMetadataCategories):
    dictOfNucCounters = {}
    for val1 in optionsForMetadataCategories[0]:
        dictOfNucCounters[val1] = NucCounter(nameOfGroup=val1, countLength=lengthOfSeqs)
        # for val2 in optionsForMetadataCategories[1]:
        #     dictOfNucCounters[val1+val2] = NucCounter(nameOfGroup=val1+val2)

            # # TODO: can have conflict if optionsForMetadataCategories[0] and optionsForMetadataCategories[1] have same values
    for val2 in optionsForMetadataCategories[1]:
        dictOfNucCounters[val2] = NucCounter(nameOfGroup=val2, countLength=lengthOfSeqs)
    return dictOfNucCounters

nucCounters = {}
# maybe better with a matrix but 2 minutes isn't that bad
# isPathogenic = np.array()
# isCow = np.array()

with open(alignedFilePath) as alignedFile:
    # for metaDataTypes in range(len(possibleClassValues)):
    #     nucCounters =
    # hard code for now


    nameOfSeq = "unknown"
    isFirstDataLine = True
    for line in alignedFile:
        line = line.strip()
        if line[0] == ">":
            nameOfSeq = line[1:]
        else:
            if isFirstDataLine:
                isFirstDataLine = False
                nucCounters = makeNucCounters(len(line), possibleClassValues)
            # nameOfNucCounter = "".join(seqNameToMetaDataType[nameOfSeq])

            for i in range(len(possibleClassValues)):
                nameOfNucCounter = seqNameToMetaDataType[nameOfSeq][i]
                charIndex = -1
                for char in line:
                    charIndex += 1
                    nucCounters[nameOfNucCounter].nucCounts[charIndex][charToIndex[char]] += 1


for metaDataCategoryIndex in range(len(possibleClassValues)):
    nucCountersOfOneMetadataCatagory = []
    for option in possibleClassValues[metaDataCategoryIndex]:
        # ["A","T","C","G"]
        # for
        nucCountersOfOneMetadataCatagory.append(nucCounters[option])
    # TODO: only works with two options in each metadatacatagory
    for nucIndex in range(len(nucCountersOfOneMetadataCatagory[0].nucCounts)):
        dist1 = nucCountersOfOneMetadataCatagory[0].getDistributionAt([nucIndex])
        dist2 = nucCountersOfOneMetadataCatagory[1].getDistributionAt([nucIndex])
        dist2 = np.array(dist2)
        dist1 = np.array(dist1)
        degOfFreedom = max(sum([val != 0 for val in dist1]), sum([val != 0 for val in dist2]))
        print(dist1, dist2, "chisquare P value??",stats.chi2.cdf(sum((dist1-dist2)**2/dist1)),df=degOfFreedom)