import re
import os
from glob import glob
# by nuc similarity (doesn't really work)
fileFolderPath = "./assemblies"
# hisC
geneToSearchSequence = "ATGAACACCGTGACTATTACCGATTTAGCGCGTGAAAACGTCCGCAACCTGACGCCGTATCAGTCAGCGCGTCGTCTGGGCGGTAACGGCGATGTCTGGCTGAACGCCAACGAATATCCCACTGCCGTGGAGTTTCAGCTCACTCAGCAAACGCTCAACCGCTACCCGGAATGCCAACCGAAAGCGGTGATTGAAAATTACGCGCAATATGCAGGCGTAAAACCGGAACAGGTGCTGGTCAGCCGTGGCGCGGACGAAGGCATTGAACTGCTGATTCGCGCTTTTTGCGAACCGGGTAAAGACGCCATCCTCTACTGCCCGCCAACGTACGGCATGTACAGCGTCAGCGCTGAAACCATTGGCGTCGAGTGCCGCACAGTGCCGACGCTGGACAACTGGCAACTGGACTTACAGGGCATTTCCGACAAGCTGGATGGCGTAAAAGTGGTCTATGTTTGCAGCCCCAACAACCCCACCGGGCAACTGATCAACCCGCAGGATTTTCGCACCCTGCTGGAGTTAACCCGCGGTAAGGCGATTGTGGTTGCCGATGAAGCTTATATCGAGTTTTGCCCGCAGGCATCGCTGGCTGGCTGGCTGGCGGAATATCCGCACCTGGCTATTTTGCGCACACTGTCGAAAGCTTTTGCTCTGGCGGGGCTTCGTTGCGGATTTACGCTGGCAAACGAAGAAGTCATCAACCTGCTGATGAAAGTGATCGCCCCCTACCCGCTCTCGACGCCGGTTGCCGACATTGCCGCACAGGCGTTAAGCCCGCAGGGAATCGTCGCCATGCGCGAACGGGTGGCGCAAATTATTGCAGAACGCGAATACCTGATTGCCGCACTGCAAGAGATCCCCTGCGTGGAGCAGGTTTTCGACTCTGAAACCAACTACATTCTGGCGCGCTTTAAAGCCTCCAGCGCAGTGTTTAAATCTTTGTGGGATCAGGGCATTATCTTACGTGATCAGAATAAACAACCCTCTTTAAGCGGCTGCCTGCGAATTACCGTCGGAACCCGTGAAGAAAGCCAGCGCGTCATTGACGCCTTACGTGCGGAGCAAGTTTGA"
geneSeq3 = "ATGAGCACCGTGACTATTACCGATTTAGCGCGTGAAAACGTCCGCAACCTGACGCCGTATCAGTCGGCGCGTCGTCTGGGCGGTAACGGCGACGTCTGGCTGAACGCCAACGAATACCCCACAGCCGTGGAGTTTCAGCTTACTCAGCAAACGCTCAACCGCTACCCGGAATGTCAGCCGAAAGCGGTGATTGAAAATTACGCGCAGTATGCAAGCGTAAAACCGGAGCAGGTGCTGGTCAGCCGTGGCGCGGACGAAGGTATTGAACTACTGATTCGCGCTTTTTGCGAACCAGGTAAAGACGCCATCCTCTACTGCCCGCCAACGTACGGCATGTACAGCGTCAGCGCTGAAACCATTGGCGTCGAGTGCCGCACAGTGCCGACGCTGGACAACTGGCAACTGGACTTGCAGGGCATTTCCGACAAGCTGGACGGCGTAAAAGTGGTCTATGTTTGCAGCCCCAACAACCCGACCGGGCAACTGATCAATCCGCAGGATTTTCGCACTCTGCTGGAGCTAACCCGCGGAAAAGCCATTGTGGTTGCCGATGAAGCCTATATTGAGTTTTGCCCGCAGGCATCGCTGGCTGGCTGGCTGGCGGAATATCCGCACCTGGCTATTTTGCGCACACTGTCGAAAGCTTTTGCTCTGGCGGGCCTTCGTTGCGGATTTACGCTGGCAAACGAAGAAGTCATCAACCTGCTGATGAAAGTGATCGCCCCCTACCCGCTCTCGACGCCGGTCGCCGACATTGCGGCCCAGGCGTTAAGCCCGCAGGGGATCGTCGCCATGCGCGAACGGGTAGCGCAAATTATTGCAGAACGCGAATATCTGATTGCCGCATTGAAAGAAATCCCCTGTGTGGAGCAGGTTTTTGACTCTGAAACCAACTACATTCTGGCGCGCTTTAAAGCCTCCAGTGCAGTGTTTAAATCTTTGTGGGATCAGGGCATTATCTTACGTGATCAGAATAAACAACCCTCTTTAAGCGGCTGCCTGCGAATTACCGTCGGAACCCGTGAAGAAAGCCAGCGCGTCATTGACGCCTTACGTGCGGAGCAAGTTTAA"

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
print(GenesAreWithinPercentIdentical(geneToSearchSequence, geneSeq3, cutoff= 0.8))
for fileName in glob(os.path.join(fileFolderPath, '*.fasta')):
    with open(fileName) as file:
        curContig = ""
        for line in file:
            if line[0] == ">":
                if len(curContig) > 50:
                    # search contig for gene
                    for i in range(len(curContig)):
                        if GenesAreWithinPercentIdentical(geneToSearchSequence, curContig[i: i+len(geneToSearchSequence)-1]):
                            print("gene found")
                curContig = ""
            else:
                curContig += line.strip()

        sameNucNumbers = []
        sameNucNumberForGene = 0
        contigIndex = 0
        reachedContigEnd = False
    print("read file")

# by gene annotations
# annotationFile1 = "AnnotationsSample1.tsv"
# annotationFile2 = "nonMastitisStrain.tsv"
# import copy
#
# matchingGenes = [] # this doesn't seem to work either
# with open(annotationFile1) as file1:
#     i = 0
#     for line1 in file1:
#         i += 1
#         cols1 = line1.split("\t")
#         # print("cols1[1]", cols1[1])
#         if cols1 == "elaB gene":
#             print("elab 1")
#         with open(annotationFile2) as file2:
#             for line2 in file2:
#                 cols2 = line2.split("\t")
#                 if cols2 == "elaB gene":
#                     print("elab 2")
#                 # print("cols2[1]",cols2[1])
#                 if cols1[1] == cols2[1] and cols1[1] != "gene":
#                     matchingGenes.append(cols1[1])
#                     # matchingGenes.append(cols2[1])
# #     print(i)
# for gene1 in copy.deepcopy(matchingGenes):
#     for _ in range(matchingGenes.count(gene1) - 1):
#             matchingGenes.remove(gene1)
# print(matchingGenes)

