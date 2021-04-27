import re
import os
from glob import glob
# by nuc similarity (doesn't really work)
refGenomeSeq = "ATGAACACCGTGACTATTACCGATTTAGCGCGTGAAAACGTCCGCAACCTGACGCCGTATCAGTCAGCGCGTCGTCTGGGCGGTAACGGCGACGTCTGGCTGAACGCCAACGAATACCCCACAGCCGTGGAGTTCCAGCTCACTCAGCAAACGCTCAACCGCTACCCGGAATGCCAGCCGAAAGCGGTGATTGAAAATTACGCGCAATATGCAGGCGTAAAACCGGAGCAGGTGCTGGTCAGCCGTGGCGCGGACGAAGGTATTGAACTGCTGATTCGCGCTTTTTGCGAACCGGGTAAAGACACCATCCTCTACTGCCCGCCAACGTACGGCATGTACAGCGTCAGCGCCGAAACCATTGGCGTCGAGTGCCGCACAGTGCCGACGCTGGACAACTGGCAACTGGACTTACAGGGCATTTCCGACAAGCTGGACGGCGTAAAAGTGGTCTATGTTTGCAGCCCCAACAACCCCACCGGGCAACTGATCAACCCGCAGGATTTTCGCACCCTGCTGGAGTTAACCCGCGGTAAGGCGATTGTGGTTGCCGATGAAGCTTATATCGAGTTTTGCCCGCAGGCATCGTTAGCTGGCTGGCTGGCGGAATATCCGCACCTGGCGATTTTGCGTACCCTTTCGAAAGCCTTCGCGTTGGCGGGGCTTCGTTGCGGATTTACGCTGGCAAACGAAGAAGTCATCAACCTGCTGATGAAAGTGATCGCCCCCTACCCGCTCTCGACGCCGGTTGCCGACATTGCCGCACAGGCGTTAAGCCCGCAGGGAATCGTCGCCATGCGCGAACGGGTGGCGCAAATTATTACAGAACGCGAATACCTGATTGCCGCACTGAAAGAGATCCCCTGCGTGGAGCAGGTTTTCGACTCTGAAACCAACTACATTCTGGCGCGCTTTAAAGCCTCCAGTGCAGTGTTTAAATCTTTGTGGGATCAGGGCATTATCTTACGTGATCAGAATAAACAACCCTCTTTAAGCGGCTGCCTGCGAATTACCGTCGGAACCCGTGAAGAAAGCCAGCGCGTCATTGACGCCTTACGTGCGGAGCAAGTTTGA"
refSeq = "ATGAACACCGTGACTATTACCGATTTAGCGCGTGAAAACGTCCGCAACCTGACGCCGTATCAGTCAGCGCGTCGTCTGGGCGGTAACGGCGATGTCTGGCTGAACGCCAACGAATATCCCACTGCCGTGGAGTTTCAGCTCACTCAGCAAACGCTCAACCGCTACCCGGAATGCCAACCGAAAGCGGTGATTGAAAATTACGCGCAATATGCAGGCGTAAAACCGGAACAGGTGCTGGTCAGCCGTGGCGCGGACGAAGGCATTGAACTGCTGATTCGCGCTTTTTGCGAACCGGGTAAAGACGCCATCCTCTACTGCCCGCCAACGTACGGCATGTACAGCGTCAGCGCTGAAACCATTGGCGTCGAGTGCCGCACAGTGCCGACGCTGGACAACTGGCAACTGGACTTACAGGGCATTTCCGACAAGCTGGATGGCGTAAAAGTGGTCTATGTTTGCAGCCCCAACAACCCCACCGGGCAACTGATCAACCCGCAGGATTTTCGCACCCTGCTGGAGTTAACCCGCGGTAAGGCGATTGTGGTTGCCGATGAAGCTTATATCGAGTTTTGCCCGCAGGCATCGCTGGCTGGCTGGCTGGCGGAATATCCGCACCTGGCTATTTTGCGCACACTGTCGAAAGCTTTTGCTCTGGCGGGGCTTCGTTGCGGATTTACGCTGGCAAACGAAGAAGTCATCAACCTGCTGATGAAAGTGATCGCCCCCTACCCGCTCTCGACGCCGGTTGCCGACATTGCCGCACAGGCGTTAAGCCCGCAGGGAATCGTCGCCATGCGCGAACGGGTGGCGCAAATTATTGCAGAACGCGAATACCTGATTGCCGCACTGCAAGAGATCCCCTGCGTGGAGCAGGTTTTCGACTCTGAAACCAACTACATTCTGGCGCGCTTTAAAGCCTCCAGCGCAGTGTTTAAATCTTTGTGGGATCAGGGCATTATCTTACGTGATCAGAATAAACAACCCTCTTTAAGCGGCTGCCTGCGAATTACCGTCGGAACCCGTGAAGAAAGCCAGCGCGTCATTGACGCCTTACGTGCGGAGCAAGTTTGA"
def GenesAreWithinPercentIdentical(seq1, seq2, cutoff= 0.8): # can do four files in <25mins - too slow
    numSame = 0
    if len(seq1) / 2 > len(seq2) or len(seq2) / 2 > len(seq1):
        return False
    for i in range(len(seq1)):
        try:
            if seq1[i] == seq2[i]:
                numSame += 1
            #if number of wrong nuc's is greater than the num possible to be within tolerance, it isn't
            if i - numSame > len(seq1)*(1-cutoff):
                return False
        except IndexError:
            print("different Gene lengths")
    return True

def translate(validSeq, codonsToAminoAcids):
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

# make codonDictionary
codonsToAminoAcids = {}
with open("codons.txt") as codonFile:
    for codonLine in codonFile:
        codonLine = codonLine.strip()
        cols = codonLine.split("\t")
        codonsToAminoAcids[cols[0]] = cols[1]

with open('ExtractionOfHisC113Genes.txt') as geneFile:
    geneSeq = ""
    numCs = 0
    numTs = 0
    for line in geneFile:
        line = line.strip()
        cols = line.split(" ")
        if re.sub(r"\d+","",cols[-1]) == "": # if only digits at the end
            # if in dna lines
            geneSeq += "".join(cols[:-1])

        elif geneSeq != "":
            print(translate(geneSeq, codonsToAminoAcids)[0:80]) # hisC mutation seems to not affect coding sequence
            if geneSeq[92] == "C":
                numCs += 1
            elif geneSeq[92] == "T":
                numTs += 1
            else:
                print("something else")
            geneSeq = ""
    if geneSeq != "": # to catch the last entry
        print(translate(geneSeq, codonsToAminoAcids)[20:40])
        if geneSeq[92] == "C":
            numCs += 1
        elif geneSeq[92] == "T":
            numTs += 1
        else:
            print("something else")
    print("%C",numCs/(numTs + numCs)) # met for
    print("%T",numTs/(numTs + numCs))
    print(numTs + numCs)