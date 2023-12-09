# from secondaryPythonScripts.functions import *
import sys

sys.path.insert(1, '/Users/cazcullimore/Documents/ericksonLabCode/secondaryPythonScripts')
try:
    from secondaryPythonScripts.functions import *
except ImportError:
    from functions import *
if len(sys.argv) < 3:
    print("please provide the paths to the gbks and the path to the file of snps sorted by significance from parsingMegacatsResults.py")
    exit(1)

pathToGBFiles = sys.argv[1]#"/Users/cazcullimore/Documents/ericksonLabCode/filesToTestMultiAlign/normalAlignedGenomes/gbks/*.gbk"
metadataFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsWithExtraMastitis.tsv"
# snpFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/snpsSortedBySignificanceWithGenesContainingThemOldIndexPathogenicity.tsv"
snpFilePath = sys.argv[2]#"/Users/cazcullimore/Documents/ericksonLabCode/RedoingEverything/snpsSortedBySignificanceWithGenesContainingThemOldIndexAnimal.tsv"
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
pathToGenes = {}

for path in glob(pathToGBFiles):
    contigs = getContigs(path)
    pathToGenes[path] = getGenesOnContigs(path, contigs)

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
        i = -1
        for path in glob(pathToGBFiles):
            i += 1
            try:
                fileName = path.split("/")[-1]
                try:
                    geneSeq = pathToGenes[path][geneName][1]
                    nuc = geneSeq[snpIndex]
                    # if i % 50 == 0:
                    print(geneSeq[snpIndex-5:snpIndex], geneSeq[snpIndex], geneSeq[snpIndex+1:snpIndex+5])
                    print(geneSeq)
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