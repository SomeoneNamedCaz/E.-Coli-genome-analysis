import re
import time
import urllib.request
import os
from glob import glob

def ReverseComplement(seq):
    reverseComplement = ''  # string of complementary nucleotides
    seq = seq.upper()  # to cover for lowercase
    for char in seq:
        if char == 'A':
            reverseComplement = 'T' + reverseComplement # add nucleotide to front (this flips the sequence)
        elif char == 'T':
            reverseComplement = "A" + reverseComplement
        elif char == 'C':
            reverseComplement = "G" + reverseComplement
        elif char == 'G':
            reverseComplement = 'C' + reverseComplement
    return reverseComplement

def GetContigs(fileName):
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

def GetGenesOnContigs(fileName, contigs):
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
                        geneSeq = ReverseComplement(geneSeq)
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

def GeneSimilarity(seq1, seq2):
    numSame = 0
    printedError = False
    for i in range(len(seq1)):
        try:
            if seq1[i] == seq2[i]:
                numSame += 1
        except IndexError:
            if not printedError:
                print("different Gene lengths")
                printedError = True
    return numSame / min(len(seq1),len(seq2))

def GenesAreWithinPercentIdentical(seq1, seq2, cutoff= 0.8):
    numSame = 0
    if len(seq1) / 2 > len(seq2) or len(seq2) / 2 > len(seq1):
        return False
    for i in range(max(len(seq1), len(seq2))):
        try:
            # if number of wrong nuc's is greater than the num possible to be within tolerance, it isn't
            if i - numSame > len(seq1)*(1-cutoff):
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
            if GenesAreWithinPercentIdentical(hypotheticalProteinSeqs[i], seqOfGeneAlreadyNamed, cutoff=cutoff):
                inList = True
        if not inList:
            nameAndSequence.append(["HP" + i, hypotheticalProteinSeqs[i]])
    return nameAndSequence
def hypotheticalProteinsAreTheSame(seq1, seq2, cutoff=0.8):
    return GenesAreWithinPercentIdentical(seq1, seq2,cutoff=0.8)

def GetSeqFilesNamesFromMasterFiles(masterFileName):
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

def DownloadSeqFileFromMasterFile(masterFileName, outputPath):
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

def StringIsOnlyATCG(string):
    string = string.upper()
    return string.count("A") + string.count("C") + string.count("T") + string.count("G") == len(string)

def FilesAreIdentical(file1Name, file2Name):
    with open(file1Name) as file1, open(file2Name) as file2:
        try:
            for line1, line2 in zip(file1,file2):
                if line1 != line2:
                    return False
        except IndexError:
            return False
    return True
def HypotheticalProteinSeqInList(seq, list, cutoff=0.8):
    for gene in list:
        if GenesAreWithinPercentIdentical(seq,gene,cutoff=cutoff):
            return True
    return False