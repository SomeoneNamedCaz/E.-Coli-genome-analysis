from secondaryPythonScripts.functions import *
from concurrent.futures import *
from time import time

if len(sys.argv) < 5:
    print("\nplease provide the path to the snp genome, the corresponding indexes, the reference sequence (.gb)")
    print("and the output path output path")
    exit(1)

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
    lastIndex = -1
    indexInSnpGenome = 0
    alignedSeq = []
    for indexOfSNPInRefSeq in indexes:
        nucsToInsert = referenceSeq[int(lastIndex) + 1:int(indexOfSNPInRefSeq)]
        alignedSeq.append(nucsToInsert)
        alignedSeq.append(seq[indexInSnpGenome])
        lastIndex = indexOfSNPInRefSeq
        indexInSnpGenome += 1

    nucsToInsert = referenceSeq[int(lastIndex) + 1:]
    alignedSeq.append(nucsToInsert)
    return "".join(alignedSeq)


numThreads = 1
t1 = time()
pool = ThreadPoolExecutor(numThreads)
futures = []
# print("starting threads")
for key in alignedGenomes.keys():
    futures.append((key,pool.submit(addRefNucs,alignedGenomes[key], snpIndexes, refSeq)))
# print("finished loading threads")

i = 0
for keyAndfuture in futures:
    alignedGenomes[keyAndfuture[0]] = keyAndfuture[1].result()
    # print("i: ", i)
    i += 1
print(numThreads,"threads took", time() - t1, "seconds")

if not os.path.exists("/".join(outputPath.split("/")[-1:])):
    os.mkdir("/".join(outputPath.split("/")[-1:]))

with open(outputPath,"w") as outFile:
    for key in alignedGenomes.keys():
        outFile.write(">" + key + "\n")
        outFile.write(alignedGenomes[key] + "\n")
        if len(refSeq) + sum([index != int(index) for index in snpIndexes]) != len(alignedGenomes[key]):
            print("PROBLEM ALIGNED GENOME IS THE WRONG SIZE")

