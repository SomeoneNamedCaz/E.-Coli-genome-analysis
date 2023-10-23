import sys
from functions import *
from collections import deque
from concurrent.futures import *
from time import time
# from dask.distributed import Client
snpGenomePath = sys.argv[1]
snpIndexesPath = sys.argv[2]
print("snp genome path", snpGenomePath)
print("index path", snpIndexesPath)

refGenomePath = sys.argv[3] # .gbff file
refGenomeContigs = getContigs(refGenomePath)
if len(refGenomeContigs) > 1:
    print("ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")

refSeq = refGenomeContigs[0]

outputPath = sys.argv[4]  # with file

snpIndexes = []
# load indexes
with open(snpIndexesPath) as file:

    for line in file:  # special stuff here to adjust for float method that was added
        line = float(line.strip())
        snpIndexes.append(line - 1)  # because genome files count starting at 1


# load snp genomes
alignedGenomes = {}
with open(snpGenomePath) as file:
    fileName = "somethingBroke" # should be written to first
    for line in file:
        line = line.strip()
        if line[0] == ">":
            fileName = line[1:]
        else:
            alignedGenomes[fileName] = line

# # TODO: make fake ref Seq for test VCFs
# snpIndexes = [0,4,5,5.1,7]
# alignedGenomes = {"genome1": deque("04510"),"genome2": deque("SEQ20")}
# refSeq = "0123456789"


def addRefNucs(seq, indexes, referenceSeq):
    refSeqDeque = deque(referenceSeq)
    # lastIndex = 0
    indexInSnpGenome = 0
    numAddedNucs = 1
    t1 = time()
    numInserts = 0
    for indexOfSNPInRefSeq in indexes:
        if indexInSnpGenome % 5000 == 0:
            print(indexInSnpGenome, time()-t1, "inserts", numInserts, "deque size", len(refSeqDeque))
            t1 = time()
            numInserts = 0
        if indexOfSNPInRefSeq - int(indexOfSNPInRefSeq) != 0:
            refSeqDeque.insert(int(indexOfSNPInRefSeq) + numAddedNucs, seq[indexInSnpGenome])
            numAddedNucs += 1
            numInserts += 1
        else:
            refSeqDeque[int(indexOfSNPInRefSeq) + numAddedNucs - 1] = seq[indexInSnpGenome]
        # nucsToInsert = refSeq[int(lastIndex):int(indexOfSNPInRefSeq)]
        # seq.insert(indexInSnpGenome, nucsToInsert)

        # lastIndex = indexOfSNPInRefSeq + 1
        indexInSnpGenome += 1

    # # add to the end of reference genome
    # nucsToInsert = refSeq[int(lastIndex):]
    # for key in alignedGenomes.keys():
    #     alignedGenomes[key].append(nucsToInsert)
    return "".join(refSeqDeque)

pool = ThreadPoolExecutor() # int(len(alignedGenomes.keys())/2)
futures = []
print("starting threads")
for key in alignedGenomes.keys():
    futures.append((key,pool.submit(addRefNucs,alignedGenomes[key], snpIndexes, refSeq)))
print("finished loading threads")
# lastIndex = 0
# indexInSnpGenome = 0
# for indexOfSNPInRefSeq in snpIndexes:
#     if indexOfSNPInRefSeq - int(indexOfSNPInRefSeq) != 0:
#         continue # if it was a deletion we can just skip that
#     nucsToInsert = refSeq[int(lastIndex):int(indexOfSNPInRefSeq)]
#     for key in alignedGenomes.keys():
#         alignedGenomes[key].insert(indexInSnpGenome, nucsToInsert)
#
#     lastIndex = indexOfSNPInRefSeq + 1
#     indexInSnpGenome += 2 # one to move to the next seq and one for the insert (since alignedGenomes is a deque of chars)

for keyAndfuture in futures:
    alignedGenomes[keyAndfuture[0]] = keyAndfuture[1].result()
    print("finished")

print(alignedGenomes)

with open(outputPath,"w") as outFile:
    for key in alignedGenomes.keys():
        outFile.write(">" + key + "\n")
        outFile.write("".join(alignedGenomes[key]) + "\n")
# print("".join(alignedGenomes["genome1"]))
