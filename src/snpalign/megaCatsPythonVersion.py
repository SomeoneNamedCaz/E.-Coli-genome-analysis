"""Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups	Results_Filename
 313 	 7.46929814283269 	 0.0238815503170304 	 2 	 Y 	 group1(19_A,_18_C,_1_T)|group2(8_A,_30_C,_1_T) 	severity-rMsaInput.txt-rResultChisqTest.txt"""
import sys
from scipy import stats
import re
import multiprocessing as mp
from functions import *
from time import time

charToIndex = {"A": 0, "T": 1, "C": 2, "G": 3, "-": 4, "N": 5}
def poolFunc(args):
    charIndex = args[0]
    alignedGenomes = args[1]
    numPossibleClassValues = args[2]
    seqNameToMetaDataType = args[3]
    nucCounters = args[4]
    # print("in func")
    for name in alignedGenomes.keys():
        for i in range(numPossibleClassValues):
            nameOfNucCounter = seqNameToMetaDataType[name][i]
            nucCounters[nameOfNucCounter][charIndex][charToIndex[alignedGenomes[name][charIndex]]] += 1
    # print("outfunc")
    
    
def calcMegaCatsStatsOnIndexWrapper(args):
    nucIndex, dists, names, matchMegacatsStyle, indexToChar, className = args
    return calcMegaCatsStatsOnIndex(nucIndex, dists, names, matchMegacatsStyle, indexToChar, className)
def calcMegaCatsStatsOnIndex(nucIndex,dists, names, matchMegacatsStyle, indexToChar,className):
    # currentNucCounter1 = nucCountersOfOneMetadataCatagory[0]
    # currentNucCounter2 = nucCountersOfOneMetadataCatagory[1]
    # dist1 = currentNucCounter1.getDistributionAt(nucIndex)
    # dist2 = currentNucCounter2.getDistributionAt(nucIndex)
    # totalGeneomesInEachMetadataGategory
    # dist2Total = sum(dist2)  # total number of strains in each catagory
    # dist1Total = sum(dist1)
    # dists = [dist1,dist2]
    distTotals = [sum(dist) for dist in dists]
    totalNumberOfGenomes = sum(distTotals)
    degOfFreedom = max(sum([val != 0 for val in dist]) for dist in dists) - 1  # num classes - 1
    if degOfFreedom == 0:  # if not a snp continue
        return ""
    if totalNumberOfGenomes == 0:
        print("zero genomes")
    chi2 = 0
    
    numsOfEachNuc = [sum([dist[i] for dist in dists]) for i in range(len(dists[0]))]  # [total number of genomes that have As, Ts,..]
    # if we have 100 genomes and 33 of the genomes have A then we expect 33*proportionOfGroup1 genomes to have A
    # print(dist1)
    # print(dist2)
    # print(nucIndex, "out of", len(nucCountersOfOneMetadataCatagory[0].nucCounts))
    distExpecteds = [[numsOfEachNuc[i] * (distTotal / totalNumberOfGenomes) for i in range(len(dists[0]))] for distTotal in distTotals]
    # dist2Expected = [dist1PlusDist2[i] * (dist2Total / totalNumberOfGenomes) for i in range(len(dist1))]
    
    for i in range(len(charToIndex.values())):
        distIndex = -1
        for distExpected in distExpecteds:
            distIndex += 1
            if distExpected[i] != 0:
                chi2 += (dists[distIndex][i] - distExpected[i]) ** 2 / distExpected[i]
        # if dist2Expected[i] != 0:
        #     chi2 += (dist2[i] - dist2Expected[i]) ** 2 / dist2Expected[i]
    
    pVal = stats.chi2.sf(chi2, df=degOfFreedom)  # for higher accuracy
    
    """Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups	Results_Filename
     313 	 7.46929814283269 	 0.0238815503170304 	 2 	 Y 	 group1(19_A,_18_C,_1_T)|group2(8_A,_30_C,_1_T) 	severity-rMsaInput.txt-rResultChisqTest.txt"""
    residueDiversity = ""
    if matchMegacatsStyle:
        indexOfDist = -1
        for dist in dists:
            megaCatsDistStr = ""
            indexOfDist += 1
            for indexWithinDist in range(len(dist)):
                if dist[indexWithinDist] != 0:
                    megaCatsDistStr += str(dist[indexWithinDist]) + "_" + indexToChar[indexWithinDist] + ","
            megaCatsDistStr = megaCatsDistStr[:-1]

            residueDiversity += names[indexOfDist] + "(" + megaCatsDistStr + ")|"
        residueDiversity = residueDiversity[:-1]
    else:
        indexOfDist = -1
        for dist in dists:
            indexOfDist += 1
            residueDiversity += names[indexOfDist] + str(dist) + "|"
        residueDiversity = residueDiversity[:-1]
    sparseTable = "N"
    if sum([sum([elt < 5 and elt != 0 for elt in dist]) > 0 for dist in dists]):
        sparseTable = "Y"
    # if pVal < 0.05:
    return "\t".join(
        [str(nucIndex + int(matchMegacatsStyle)), str(chi2), str(pVal), str(degOfFreedom), sparseTable,
         residueDiversity, className]) + "\n"
    # return ""

def calculateMegaCatsStats(alignedFilePath, metadataFilePath, outFilePath, numProcesses=16, matchMegacatsStyle=True, debug=False):
    
    indexToChar = {}
    for k,v in charToIndex.items():
        indexToChar[v] = k
    
    
    classNames = []
    possibleClassValues = [] # maybe do a separate class for this
    seqNameToMetaDataType = {}
    with open(metadataFilePath) as metadataFile: # process metadata
        isFirstLine = True
        for line in metadataFile:
            line = line.strip()
            cols = line.split("\t")
            if isFirstLine:
                classNames = cols[1:]
                for _ in range(len(cols)-1):
                    possibleClassValues.append(set())
                isFirstLine = False
            else:
                metadataType = cols[1:]
                for colIndex in range(len(cols[1:])):
                    possibleClassValues[colIndex].add(cols[colIndex + 1])
                    metadataType[colIndex] += classNames[colIndex]
                seqNameToMetaDataType[re.sub("\..+", "", cols[0])] = metadataType
    
    if debug:
        print(possibleClassValues)
    
    nucCounters = {}
    alignedGenomes = {}
    argsForCountNucs = []
    lineLen = -1
    t1 = time()
    
    nucCounters = countNucs(readInFastaAsList(alignedFilePath),seqNameToMetaDataType,possibleClassValues,classNames)
    if debug:
        print("took",time()-t1)
    
    with open(outFilePath, "w") as outFile:
        headerInfo = "Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups	metaDataCategory"
        outFile.write(headerInfo + "\n")
        for metaDataCategoryIndex in range(len(possibleClassValues) - 1, -1, -1):
            
            nucCountersOfOneMetadataCatagory = []
            metaDataCategoryName = classNames[metaDataCategoryIndex]
            try:
                for option in possibleClassValues[metaDataCategoryIndex]:
                    
                    # ["A","T","C","G"]
                    nucCountersOfOneMetadataCatagory.append(nucCounters[option + metaDataCategoryName])
            except KeyError:
                print("WARNING: couldn't find genomes with metadata for a category")
                continue
            # TODO: only works with two options in each metadatacatagory
            
            args = []
            for nucIndex in range(len(nucCountersOfOneMetadataCatagory[0].nucCounts)):
                dists = []
                names = []
                for counter in nucCountersOfOneMetadataCatagory:
                    dists.append(counter.getDistributionAt(nucIndex))
                    names.append(counter.nameOfGroup)
                args.append((nucIndex,dists,names,matchMegacatsStyle,indexToChar,classNames[metaDataCategoryIndex]))
            pool = mp.Pool(numProcesses)
            stuffToWrite = pool.map(calcMegaCatsStatsOnIndexWrapper, args)
            for elt in stuffToWrite:
                outFile.write(elt)
    
class NucCounter():
    """
        nuc counter is a class to hold the distribution of nucleotides of all genomes of a certain metadata category at each position
        where each element of the outer list represents a different position in the genome
        and each inner list contains the number of genomes that have a particular nulcleotide at that position
        i.e. [[num_genomes that have A at position 1, num_genomes that have T at position 1, num Cs, num Gs, num "-"s, num Ns], [As at position 2 ...] ...]

    """
    
    def __init__(self, countLength=0, nameOfGroup="unnamed", ):
        self.nucCounts = []
        for _ in range(countLength):
            self.nucCounts.append([0] * len(charToIndex.values()))
        self.nameOfGroup = nameOfGroup
    
    def getDistributionAt(self, positionInGenome):
        return self.nucCounts[positionInGenome]
    
    def appendDistribution(self, distribution):
        self.nucCounts.append(distribution)
        
    def __repr__(self):
        return  '(' + str(self.nucCounts) + ', ' + str(self.nameOfGroup) + ')'
    def __str__(self):
        return self.__repr__()
def makeNucCounters(lengthOfSeqs, optionsForMetadataCategories, metadataCategoryTitles):
    dictOfNucCounters = {}
    for metadataCategoryIndex in range(len(metadataCategoryTitles)):

        for val1 in optionsForMetadataCategories[metadataCategoryIndex]:
            # this adds the category title so that you can have redundant naming schemes for each metadata category and not have conflicts
            dictOfNucCounters[val1 + metadataCategoryTitles[metadataCategoryIndex]] = NucCounter(nameOfGroup=val1, countLength=lengthOfSeqs)
        
        # # TODO: can have conflict if optionsForMetadataCategories[0] and optionsForMetadataCategories[1] have same values
    # for val2 in optionsForMetadataCategories[1]:
    #     dictOfNucCounters[val2] = NucCounter(nameOfGroup=val2, countLength=lengthOfSeqs)
    return dictOfNucCounters
def countNucs(alignedFileList,seqNameToMetaDataType,possibleClassValues,metadataCategoryTitles, debug=False):
    nameOfSeq = "unknown"
    isFirstDataLine = True
    for line in alignedFileList:
        line = line.strip()
        if line == "":
            continue
        if line[0] == ">":
            
            if debug:
                print("before", line[1:])
            nameOfSeq = re.sub("\..+", "", line[1:])
            
            if debug:
                print("nameOfSeq", nameOfSeq)
        else:
            if nameOfSeq == nameOfRefSnpGenome:
                continue
            if isFirstDataLine:
                isFirstDataLine = False
                nucCounters = makeNucCounters(len(line), possibleClassValues,metadataCategoryTitles)
            
            for i in range(len(possibleClassValues)):
                if nameOfSeq == nameOfRefSnpGenome:
                    continue
                nameOfNucCounter = seqNameToMetaDataType[nameOfSeq][i]
                charIndex = -1
                for char in line:
                    charIndex += 1
                    nucCounters[nameOfNucCounter].nucCounts[charIndex][charToIndex[char]] += 1
    return nucCounters
                
if __name__ == "__main__":
    snpAlignmentFilePath = sys.argv[1]
    metadataFilePath = sys.argv[2]
    outFilePath = sys.argv[3]
    matchMegacatsStyle = bool(sys.argv[4])
    calculateMegaCatsStats(snpAlignmentFilePath, metadataFilePath, outFilePath, matchMegacatsStyle)