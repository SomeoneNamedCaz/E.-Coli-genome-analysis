import re
import time
import urllib.request
import os
from glob import glob
import threading
import sys

def reverseComplement(seq):
    reverseComplementSeq = ''  # string of complementary nucleotides
    seq = seq.upper()  # to cover for lowercase
    for char in seq:
        if char == 'A':
            reverseComplementSeq = 'T' + reverseComplementSeq # add nucleotide to front (this flips the sequence)
        elif char == 'T':
            reverseComplementSeq = "A" + reverseComplementSeq
        elif char == 'C':
            reverseComplementSeq = "G" + reverseComplementSeq
        elif char == 'G':
            reverseComplementSeq = 'C' + reverseComplementSeq
    return reverseComplementSeq

def getContigs(fileName):
    contigs = []
    with open(fileName) as file:
        contig = ""
        inContig = False
        for line in file:
            line = line.strip()
            cols = line.split()
            if line == 'ORIGIN':
                inContig = True
            elif line == "//":
                contig = contig.upper()
                contigs.append(contig)
                contig = ""
                inContig = False
            if inContig:
                contig += "".join(cols[1:]) # everything but numbers
    return contigs

def getGenesOnContigs(fileName, contigs): #TODO: doesn't do such a good job on <,>, or num..num,num..num
    genes = {} # dictionary keys gene_name values = products and seqs
    with open(fileName) as file:
        contig = ""
        inContig = False
        inFeatures = False
        inCDS = False
        geneStartPoses = []
        geneStopPoses = []
        geneIsForward = []
        geneNames = []
        geneProducts = []
        geneSeq = ""
        contigIndex = -1
        numHypotheticalProteins = 0
        for line in file:
            line = line.strip()
            cols = line.split()
            if len(cols) == 0:
                continue
            if cols[0] == "FEATURES":
                inFeatures = True
            elif line == 'ORIGIN':
                contigIndex += 1
                inContig = True
                inCDS = False
                # if geneStartPoses == [] or geneStopPoses == []:
                #     raise Exception("ERROR: found contig before annotations")
            elif line == "//":
                contig = ""
                inContig = False
            if inContig: # if in contig
                # print("Startpos",len(geneStartPoses))
                for geneIndex in range(len(geneStartPoses)):  # for each gene
                    # print("index", geneIndex)
                    geneStartPos = geneStartPoses[geneIndex]
                    geneStopPos = geneStopPoses[geneIndex]
                    geneSeq = contigs[contigIndex][geneStartPos - 1:geneStopPos]
                    geneSeq = geneSeq.upper()
                    if not geneIsForward[geneIndex]:
                        geneSeq = reverseComplement(geneSeq)
                    try:
                        if geneNames[geneIndex] != "unnamed": # use names as keys if you can
                            genes[geneNames[geneIndex]] = [geneProducts[geneIndex], geneSeq]
                        else:
                            if geneProducts[geneIndex] == "hypothetical protein":
                                genes[geneSeq] = [geneProducts[geneIndex], geneSeq]
                                numHypotheticalProteins += 1
                            else:
                                genes[geneProducts[geneIndex]] = [geneProducts[geneIndex], geneSeq]
                    except IndexError:
                        0#print("noname or product")
                geneStartPoses = []
                geneStopPoses = []
                geneIsForward = []
                geneNames = []
                geneProducts = []
            elif inFeatures: # if in annotations
                if cols[0] == "CDS": # when in gene line
                    geneNames.append("unnamed")
                    inCDS = True
                    nums = re.sub(r"complement", "", cols[1])
                    if nums != cols[1]:
                        geneIsForward.append(False)
                    else:
                        geneIsForward.append(True)
                    shortenEndOfGenomeBy = 0
                    if ">" in nums:
                        shortenEndOfGenomeBy = 1
                    nums = re.sub(r"[\(\)<>]", "", nums) # remove perentheses
                    nums = nums.split("..")
                    geneStartPoses.append(int(nums[0]))
                    geneStopPoses.append(int(nums[1])-shortenEndOfGenomeBy)
                elif cols[0] == "gene":
                    inCDS = False
                elif inCDS:
                    if line[:7] == '/gene="':
                        geneNames[-1] = re.sub('/gene="', "", line)[:-1]
                    elif line[:10] == '/product="':
                        geneProducts.append(re.sub('/product="', "", line)[:-1])
    return genes
class SNP:
    def __init__(self, location, oldNuc, newNuc, pValue):
        self.location = location
        self.oldNuc = oldNuc
        self.newNuc = newNuc
        self.pValue = pValue

class Gene:
    def __init__(self, startPos, stopPos, sequence, name, product):
        self.startPos = startPos
        self.stopPos = stopPos
        self.sequence = sequence
        self.name = name
        self.product = product
        self.snps = [] # list of SNP objects
        self.counter = 0


def getGenesOnContigsByPosition(fileName, contigs):
    """":returns list of gene info ordered by position elements are of the gene class"""
    genes = []
    with open(fileName) as file:
        contig = ""
        inContig = False
        inFeatures = False
        inCDS = False
        geneStartPoses = []
        geneStopPoses = []
        geneIsForward = []
        geneNames = []
        geneProducts = []
        geneSeq = ""
        contigIndex = -1
        numHypotheticalProteins = 0
        for line in file:
            line = line.strip()
            cols = line.split()
            if len(cols) == 0:
                continue
            if cols[0] == "FEATURES":
                inFeatures = True
            elif line == 'ORIGIN':
                contigIndex += 1
                inContig = True
                inCDS = False
                # if geneStartPoses == [] or geneStopPoses == []:
                #     raise Exception("ERROR: found contig before annotations")
            elif line == "//":
                contig = ""
                inContig = False
            if inContig: # if in contig
                # print("Startpos",len(geneStartPoses))
                for geneIndex in range(len(geneStartPoses)):  # for each gene
                    # print("index", geneIndex)
                    geneStartPos = geneStartPoses[geneIndex]
                    geneStopPos = geneStopPoses[geneIndex]
                    geneSeq = contigs[contigIndex][geneStartPos - 1:geneStopPos]
                    geneSeq = geneSeq.upper()
                    if not geneIsForward[geneIndex]:
                        geneSeq = reverseComplement(geneSeq)
                    try:
                        genes.append(Gene(geneStartPos, geneStopPos, geneSeq, geneNames[geneIndex], geneProducts[geneIndex]))
                        if geneProducts[geneIndex] == "hypothetical protein":
                            numHypotheticalProteins += 1
                    except IndexError:

                        print("noname or product")
                geneStartPoses = []
                geneStopPoses = []
                geneIsForward = []
                geneNames = []
                geneProducts = []
            elif inFeatures: # if in annotations
                if cols[0] == "CDS": # when in gene line
                    geneNames.append("unnamed")
                    geneProducts.append("no product")
                    inCDS = True
                    nums = re.sub(r"complement", "", cols[1])
                    if nums != cols[1]:
                        geneIsForward.append(False)
                    else:
                        geneIsForward.append(True)
                    shortenEndOfGenomeBy = 0
                    if ">" in nums:
                        shortenEndOfGenomeBy = 1
                    # TODO: include if the gene is a partial CDS (has a < or >)
                    nums = re.sub(r"[A-Za-z\(\)<>]+", "", nums) # remove perentheses
                    # if "," in nums: # if of form join(this..this, this..this)

                    nums = nums.split("..")
                    geneStartPoses.append(int(nums[0]))
                    geneStopPoses.append(int(nums[-1])-shortenEndOfGenomeBy)
                elif cols[0] == "gene":
                    inCDS = False
                elif inCDS:
                    if line[:7] == '/gene="':
                        geneNames[-1] = re.sub('/gene="', "", line)[:-1]
                    elif line[:10] == '/product="':
                        geneProducts[-1] = re.sub('/product="', "", line)[:-1]
    genes.sort(key=lambda a: a.startPos)
    return genes

def geneSimilarity(seq1, seq2):
    numSame = 0
    printedError = False
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    for i in range(len(seq1)):

        try:
            if seq1[i] == seq2[i]:
                numSame += 1
        except IndexError:
            if not printedError:
                print("different Gene lengths")
                printedError = True
    return numSame / max(len(seq1),len(seq2))

def genesAreWithinPercentIdentical(seq1, seq2, cutoff= 0.8):
    numSame = 0
    if len(seq1) / 2 > len(seq2) or len(seq2) / 2 > len(seq1):
        return False
    for i in range(max(len(seq1), len(seq2))):
        try:
            # if number of wrong nuc's is greater than the num possible to be within tolerance, it isn't
            if i - numSame > max(len(seq1), len(seq2))*(1-cutoff):
                return False
            if seq1[i] == seq2[i]:
                numSame += 1
        except IndexError:
            0
    return True

def SNPLocations(seq1, seq2):
    locationsOfDifferences = []
    for i in range(len(seq1)):
        try:
            if seq1[i] != seq2[i]:
                locationsOfDifferences.append(i)
        except IndexError:
            print("different Gene lengths")
    return locationsOfDifferences

def translate(validSeq):
    codonsToAminoAcids = makeCodonDictionary()
    RNASeq = re.sub("T", "U", validSeq)
    currCodon = ''
    aminoAcidSeq = ''
    for char in RNASeq:
        currCodon += char
        if len(currCodon) == 3:
            aminoAcid = codonsToAminoAcids[currCodon]
            currCodon = ''
            if aminoAcid == '*':
                break
            else:
                aminoAcidSeq += aminoAcid
    return aminoAcidSeq

def makeCodonDictionary():
    codonsToAminoAcids = {}
    with open("codons.txt") as codonFile:
        for codonLine in codonFile:
            codonLine = codonLine.strip()
            cols = codonLine.split("\t")
            codonsToAminoAcids[cols[0]] = cols[1]
    return codonsToAminoAcids

def convertHypotheticalProteinsIntoTempNames(hypotheticalProteinSeqs, cutoff=0.8):
    nameAndSequence = []
    for i in range(len(hypotheticalProteinSeqs)):
        inList = False
        for geneAlreadyNamed in nameAndSequence:
            seqOfGeneAlreadyNamed = geneAlreadyNamed[1]
            if genesAreWithinPercentIdentical(hypotheticalProteinSeqs[i], seqOfGeneAlreadyNamed, cutoff=cutoff):
                inList = True
        if not inList:
            nameAndSequence.append(["HP" + str(i), hypotheticalProteinSeqs[i]])
    return nameAndSequence
def hypotheticalProteinsAreTheSame(seq1, seq2, cutoff=0.8):
    return genesAreWithinPercentIdentical(seq1, seq2, cutoff=0.8)

def getSeqFilesNamesFromMasterFiles(masterFileName):
    """
    :param masterFileName:
    :return: list of seqFileNames that the master file refers to
    """
    seqFilesNames = []
    onFirstLine = True
    isMasterRecord = False
    with open(masterFileName) as masterFile:
        for line in masterFile:
            line = line.strip()
            cols = line.split()
            if onFirstLine and cols[1][-1] == "0":  # if master record
                isMasterRecord = True
            if isMasterRecord and cols[0] == "WGS":
                startFileName = cols[1].split("-")[0]
                stopFileName = cols[1].split("-")[1]
                numberOfFiles = int(stopFileName[-5:])
                for i in range(numberOfFiles):
                    indexAsString = str(i + 1)
                    # add to 5 len
                    while len(indexAsString) < 5:
                        indexAsString = "0" + indexAsString
                    contigName = startFileName[:-5] + indexAsString
                    seqFilesNames.append(contigName)
                    # t1 = time.time()
                    # url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + contigName + "&retType=gb"
                    # urllib.request.urlretrieve(url, "./DownloadedFromSSB/" + contigName + ".gb")
                    # time.sleep(max(1 / 3 + extraWaitTime - (time.time() - t1), 0))
            onFirstLine = False
        return seqFilesNames

def downloadSeqFileFromMasterFile(masterFileName, outputPath):
    seqFilesNames = []
    onFirstLine = True
    isMasterRecord = False
    with open(masterFileName) as masterFile:
        for line in masterFile:
            cols = line.split()
            line = line.strip()
            if len(cols) == 0:
                continue
            if onFirstLine and cols[1][-1] == "0":  # if master record
                isMasterRecord = True
            if isMasterRecord and cols[0] == "WGS":
                startFileName = cols[1].split("-")[0]
                stopFileName = cols[1].split("-")[1]
                numberOfFiles = int(stopFileName[-5:])
                for i in range(numberOfFiles):
                    indexAsString = str(i + 1)
                    # add to 5 len
                    while len(indexAsString) < 5:
                        indexAsString = "0" + indexAsString
                    contigName = startFileName[:-5] + indexAsString
                    t1 = time.time()
                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + contigName + "&retType=gb"
                    try:
                        urllib.request.urlretrieve(url, outputPath + contigName + ".gb")
                    except:
                        time.sleep(2.2/3)
                        urllib.request.urlretrieve(url, outputPath + contigName + ".gb")

                    time.sleep(max(2 / 3 - (time.time() - t1), 0))
            onFirstLine = False

def stringIsOnlyATCG(string):
    string = string.upper()
    return string.count("A") + string.count("C") + string.count("T") + string.count("G") == len(string)

def filesAreIdentical(file1Name, file2Name):
    with open(file1Name) as file1, open(file2Name) as file2:
        try:
            for line1, line2 in zip(file1,file2):
                if line1 != line2:
                    return False
        except IndexError:
            return False
    return True
def hypotheticalProteinSeqInList(seq, list, cutoff=0.8):
    for gene in list:
        if genesAreWithinPercentIdentical(seq, gene, cutoff=cutoff):
            return True
    return False

def compareGenesByAAGroup(seq1, seq2, cutoff=0.8):
    protein1 = translate(seq1)
    protein2 = translate(seq2)
    #          pos polar        neg polar   polar                   nonpolar
    aaGroups = [["R","H","K"],["D","E"], ["S", "T", "N", "Q"], ["A", "V", "I", "L", "M", "F", "Y", "W"]]
    numSame = 0
    for aaIndex in range(len(protein1)):
        aa1 = protein1
        aa2 = protein2
        if aa1 == aa2:
            numSame += 1
        else:
            for aaGroup in aaGroups:
                if aa1 in aaGroup and aa2 in aaGroup:
                    numSame += 1
                    break
    precentSame =  numSame / max(len(protein1), len(protein2))
    if precentSame >= cutoff:
        return True
    return False

def filesContainTheSameInformation(fileName1, fileName2):
    with open(fileName1) as f1:
        for l1 in f1:
            foundLine = False
            with open(fileName2) as f2:
                for l2 in f2:
                    if l1 == l2:
                        foundLine = True
            if not foundLine:
                return False
    return True

def separateGenomesFromSingleFastaFile(filePath):
    # the file needs to have >LOCUS_NAME other stuff for each contig
    # LOCUS_Name needs to be one word
    # sequences can span multiple lines
    # other stuff needs to contain "complete genome" if it is a complete genome or "whole genome shotgun sequence"
    pathLocation = '/'.join(filePath.split("/")[:-1])
    with open(filePath) as fileToSort:
        onFirstLine = True
        lastLocusName = ""
        lastStrainName = ""
        lastContigSeq = ""
        lastHeaderLine = ""
        lastStrainName = ""
        for line in fileToSort:
            if line.strip() == "":
                continue
            if line[0] == ">":
                if not onFirstLine:
                    genomeFilePath = pathLocation + "/" + lastStrainName + ".fasta"
                    fileData = []
                    try:
                        with open(genomeFilePath) as fileToWrite:
                            for data in fileToWrite:
                                fileData.append(data)
                    except FileNotFoundError:
                        0  # print("make newFile")
                    with open(genomeFilePath, "w") as fileToWrite:
                        for data in fileData:
                            fileToWrite.write(data)
                        fileToWrite.write(lastHeaderLine)
                        fileToWrite.write(lastContigSeq)
                        lastContigSeq = ""
                        lastHeaderLine = ""
                lastLocusName = line[1:].split()[0]  # first word without ">"
                if "complete genome" in line:
                    lastStrainName = lastLocusName
                elif "whole genome shotgun sequence" in line:
                    lastStrainName = lastLocusName[:-5]
                else:
                    print("error not shotgun or complete")
                    lastStrainName = "weird_" + lastLocusName[:-5]

                lastHeaderLine = line
                onFirstLine = False
            else:
                lastContigSeq += line  # no strip to keep formatting

def aGeneIsAPartofAnotherGene(gene1, gene2, cutoff=0.8):
    if len(gene2) > len(gene1):
        largerGene = gene2
        smallerGene = gene1
    else:
        largerGene = gene1
        smallerGene = gene2

    for startPosition in range(len(largerGene) - len(smallerGene)):
        numSame = 0
        for nucIndex in range(len(smallerGene)):
            if smallerGene[nucIndex] == largerGene[startPosition + nucIndex]:
                numSame += 1
        if numSame/len(smallerGene) > cutoff:
            return True
    return False

def aGeneIsAPartofAnotherGeneFast(gene1, gene2, cutoff=0.8):
    if len(gene2) > len(gene1):
        largerGene = gene2
        smallerGene = gene1
    else:
        largerGene = gene1
        smallerGene = gene2
    listOfComparisons = []
    returnFlag = False
    class compareGeneThread(threading.Thread):
        def __init__(self, genePart1, genePart2):
            threading.Thread.__init__(self)
            self.genePart1 = genePart1
            self.genePart2 = genePart2
        def run(self):
            if (genesAreWithinPercentIdentical(self.genePart1, self.genePart2, cutoff=cutoff)):
                listOfComparisons.append(True)
    for startPosition in range(len(largerGene) - len(smallerGene)):
        if startPosition % 100 == 0:
            print(startPosition)
        compareGeneThread(largerGene[startPosition:], smallerGene).start()
        # numSame = 0
        # if (genesAreWithinPercentIdentical(largerGene[startPosition:], smallerGene, cutoff=cutoff)):
        #     return True

    # geneThread.join()
    # time.sleep(10) # for if one thread isn't finished later
    # for comparison in listOfComparisons:
    #     if comparison:
    #         return True
    return returnFlag

def aGeneIsAPartofAnotherGeneFastST(gene1, gene2, cutoff=0.8):
    if len(gene2) > len(gene1):
        largerGene = gene2
        smallerGene = gene1
    else:
        largerGene = gene1
        smallerGene = gene2
    listOfComparisons = []
    for startPosition in range(len(largerGene) - len(smallerGene)):
        if startPosition % 100 == 0:
            print(startPosition)
        if (genesAreWithinPercentIdentical(largerGene[startPosition:], smallerGene, cutoff=cutoff)):
            return True
    return False
def insertAtIndex(stringToInsert, index, char):
    """ returns the string that has the inserted index"""
    stringToInsert = stringToInsert[:index] + char + stringToInsert[index:]
    return stringToInsert