from secondaryPythonScripts.functions import *
from concurrent.futures import *
from time import time


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
def reconstructNormalAlignment(snpGenomePath, snpIndexesPath, refGenomePath, outputPath):
    print("snp genome path", snpGenomePath)
    print("index path", snpIndexesPath)
    print("refGenomePath", refGenomePath)
    print("output path", outputPath)
    refGenomeContigs = getContigs(refGenomePath)
    if len(refGenomeContigs) > 1:
        print(
            "ref genome has more than one contig, this is not a problem as long as the first contig is the chromosome")
    refSeq = refGenomeContigs[0]
    
    snpIndexes = []
    # load indexes
    with open(snpIndexesPath) as file:
    
        for line in file:  # special stuff here to adjust for float method that was added
            line = float(line.strip())
            snpIndexes.append(line)
    
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

    
    numThreads = 1
    t1 = time()
    pool = ThreadPoolExecutor(numThreads)
    futures = []
    
    for key in alignedGenomes.keys():
        futures.append((key,pool.submit(addRefNucs,alignedGenomes[key], snpIndexes, refSeq)))
    
    i = 0
    for keyAndfuture in futures:
        alignedGenomes[keyAndfuture[0]] = keyAndfuture[1].result()
        # print("i: ", i)
        i += 1
    print(numThreads,"threads took", time() - t1, "seconds")
    
    if not os.path.exists("/".join(outputPath.split("/")[-1:])) and outputPath.count("/") != 0:
        os.mkdir("/".join(outputPath.split("/")[-1:]))
    
    with open(outputPath,"w") as outFile:
        for key in alignedGenomes.keys():
            outFile.write(">" + key + "\n")
            outFile.write(alignedGenomes[key] + "\n")
            if len(refSeq) + sum([index != int(index) for index in snpIndexes]) != len(alignedGenomes[key]):
                print("PROBLEM ALIGNED GENOME IS THE WRONG SIZE")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("\nplease provide the path to the snp genome, the corresponding indexes, the reference sequence (.gb)")
        print("and the output path output path")
        exit(1)
    
    snpGenomePath = sys.argv[1]
    snpIndexesPath = sys.argv[2]
    print("snp genome path", snpGenomePath)
    print("index path", snpIndexesPath)
    
    refGenomePath = sys.argv[3]  # .gbff file
    
    outputPath = sys.argv[4]  # with file
    reconstructNormalAlignment(snpGenomePath, snpIndexesPath, refGenomePath, outputPath)