import sys
import time
from concurrent.futures import *
from glob import glob

sys.path.insert(1, '/scr/snpalign/secondaryPythonScripts')
try:
    from scr.snpalign.functions import *
except ImportError:
    pass
from functions import *
if len(sys.argv) < 3:
    print("please provide the paths to the gbks and the path to the phylogroupSnpFile of snps sorted by significance from parsingMegacatsResults.py")
    exit(1)

pathToGBFiles = sys.argv[1]#DATA_DIR + "filesToTestMultiAlign/normalAlignedGenomes/gbks/*.gbk"
metadataFilePath = DATA_DIR + "metaDataForMetaCatsWithExtraMastitis.tsv"
snpFilePath = sys.argv[2]#DATA_DIR + "RedoingEverything/snpsSortedBySignificanceWithGenesContainingThemOldIndexAnimal.tsv"
annotatedRefGenomePath = DATA_DIR + "refGenomes/k-12.gbff"
print("gb files",pathToGBFiles)
print("metadata", metadataFilePath)
print("megaCatsFile",snpFilePath)
fileNameToMetadataCategory = {}
with open(metadataFilePath) as metadataFile: # process metadata
    isFirstLine = True
    for line in metadataFile:
        line = line.strip()
        cols = line.split("\t")
        if isFirstLine:
            isFirstLine = False
        else:
            fileNameToMetadataCategory[cols[0].split(".")[0]] = cols[1:]
# outFile = "qorB."
# snpIndex = 230
# # for 56 {'pathogen': {'A': 0, 'T': 148, 'C': 0, 'G': 335}, 'commensal': {'A': 3, 'T': 226, 'C': 0, 'G': 233}, 'cow': {'A': 0, 'T': 255, 'C': 0, 'G': 141}, 'chicken': {'A': 3, 'T': 119, 'C': 0, 'G': 427}}
# geneName = "thrB"
# print(geneName, snpIndex)
numThreads = 64
gbFiles = glob(pathToGBFiles)[::10][:10] # to get evenly spaced genomes

pool = ThreadPoolExecutor(min(numThreads, len(gbFiles)))
futures = []
pathToGenes = {}

for path in gbFiles:
    pathToGenes[path] = {}
t1 = time.time()

def x(path):
    contigs = getContigs(path)
    return getGenesOnContigs(path, contigs)
for path in gbFiles:
    futures.append((path,pool.submit(x, path)))
    
print("started")
refGenes = getGenesOnContigs(annotatedRefGenomePath, getContigs(annotatedRefGenomePath))
i = 0
for key, future in futures:
    i += 1
    pathToGenes[key] = future.result()
    print(i, "done")
print("finished")
print(numThreads, " threads took ", time.time()-t1,"seconds")

with open(snpFilePath) as snpFile:
    isFirstLine = True
    for line in snpFile:
        if isFirstLine:
            isFirstLine = False
            continue
        geneDistributions = {"pathogen": {"A":0, "T":0, "C":0, "G":0}, "commensal": {"A":0, "T":0, "C":0, "G":0}, "cow": {"A":0, "T":0, "C":0, "G":0}, "chicken": {"A":0, "T":0, "C":0, "G":0}}
        cols = line.split("\t")
        geneName = cols[3]
        snpIndex = int(cols[2])
        try:
            refSeq = refGenes[geneName].sequence
            if not refGenes[geneName].isForward:
                refSeq = reverseComplement(refSeq)
        except KeyError:
            continue
        i = -1
        for path in gbFiles:
            i += 1
            try:
                fileName = path.split("/")[-1]
                try:
                    gene = pathToGenes[path][geneName]
                    geneSeq = gene.sequence
                    if not gene.isForward:
                        geneSeq = reverseComplement(geneSeq)
                    nuc = geneSeq[snpIndex]
                    # if i % 50 == 0:
                    print("scaf forward:",end="")
                    print(gene.isForward, "\t", geneSeq[snpIndex - 5:snpIndex],
                          geneSeq[snpIndex],
                          geneSeq[snpIndex + 1:snpIndex + 5])
                    print("ref  forward:", end="")
                    print(refGenes[geneName].isForward, "\t", refSeq[snpIndex - 5:snpIndex],
                          refSeq[snpIndex],
                          refSeq[snpIndex + 1:snpIndex + 5])
                    print("same length?",len(refGenes[geneName].sequence) == len(gene.sequence))
                    secondMetadataCatagoryForFile = fileNameToMetadataCategory["scaffold_" + fileName.split(".")[0]][1]
                    firstMetadataCatagoryForFile = fileNameToMetadataCategory["scaffold_" + fileName.split(".")[0]][0]
                    try:
                        geneDistributions[secondMetadataCatagoryForFile][nuc] += 1
                        geneDistributions[firstMetadataCatagoryForFile][nuc] += 1
                    except:
                        print(nuc, "not in metadata")
                except KeyError:
                    print(fileName, "didn't have the gene", geneName)
            except KeyboardInterrupt:
                exit(0)
            except Exception as e:
                print(e)
        print(geneDistributions)
        print(line)
