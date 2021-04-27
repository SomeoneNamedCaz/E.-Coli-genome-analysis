import sys
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
    return numSame / min(len(seq1), len(seq2))

geneSeqFileName = sys.argv[1]
geneSeqs = []
with open(geneSeqFileName) as geneSeqFile:
    for line in geneSeqFile: # each line is a new gene seq in this file
        geneStart = line.find(":")+1
        gene = line[geneStart:]
        geneSeqs.append(gene.strip())

for index1 in range(len(geneSeqs)):
    for index2 in range(index1 + 1, len(geneSeqs)): # so no repeat comparisions
        print(index1, "v.", index2,":", GeneSimilarity(geneSeqs[index1], geneSeqs[index2]))

# codon = []
# idx = 0
# for nuc in geneSeq3:
#     codon.append(nuc)
#     if len(codon) == 3:
#         print(idx, codon)
#         codon = []
#     idx += 1
