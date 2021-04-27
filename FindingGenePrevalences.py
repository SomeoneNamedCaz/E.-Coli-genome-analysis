import os
from glob import glob
import re
from functions import *
import copy

outFileName = "genePrevalencesMastitis.tsv" # structured as geneName + \t + count
pathToGenomes = "./annotatedGenomes/"
# pathToGenomes = "./natureStuff/annotatedPathogenNatureGenomes/gbks/"
# pathToGenomes = "./natureStuff/annotatedCommensalNatureGenomes/gbks/"
# pathToGenomes = "./Annotations/DownloadingFilesFromNCBI/AllCommensalStrains/assemblies/bovineCommensalAnnotatedGenomes/gbks/"
writtenToOutFileBefore = False
fileCount = 0
for annotationFileName in glob(pathToGenomes + "*"):
    fileCount += 1
    print(fileCount)
    contigs = GetContigs(annotationFileName)
    genes = GetGenesOnContigs(annotationFileName, contigs)
    fileData = []
    if writtenToOutFileBefore: # read in data
        with open(outFileName) as outFile:
            for line in outFile:
                fileData.append(line)
                # geneName = line.split("\t")[0]
    # write out new counts
    with open(outFileName, "w") as outFile:
        for line in fileData:
            cols = line.split("\t")
            geneName = cols[0]
            geneCount = int(cols[1])
            # print(geneCount)
            if geneName in genes.keys():
                # print("added")
                geneCount += 1
                del genes[geneName]

            # if looking at hypothetical protein line in the count file
            elif geneName.count("A") + geneName.count("C") + geneName.count("T") + geneName.count("G") == len(geneName):
                foundSeq = ""
                for hypotheticalProteinSeq in genes.keys():
                    if GenesAreWithinPercentIdentical(hypotheticalProteinSeq, geneName, cutoff=0.8):
                        geneCount += 1
                        foundSeq = hypotheticalProteinSeq
                if foundSeq != "":
                    del genes[foundSeq]
            outFile.write(geneName + "\t" + str(geneCount) + "\n")
        # write out new genes present
        for geneName in genes.keys():
            outFile.write(geneName + "\t1\n")
        writtenToOutFileBefore = True







    # with open(annotationFileName) as annotationFile:
    #     for line in annotationFile:

# withlist
# geneCounts = [] # of [name, count]
# for annotationFileName in glob(pathToGenomes + "*"):
#     fileCount += 1
#     print(fileCount)
#     contigs = GetContigs(annotationFileName)
#     genes = GetGenesOnContigs(annotationFileName, contigs)
#     fileData = []
#     for geneToAdd in genes:
#         addedToExistingGenes = False
#
#         for geneNameAndCount in geneCounts:
#             name = geneNameAndCount[0]
#             count = geneNameAndCount[1]
#             if name == geneToAdd or GenesAreWithinPercentIdentical(name, geneToAdd):
#                 count += 1
#                 addedToExistingGenes = True
#                 # print("wrong")
#         if not addedToExistingGenes:
#             geneCounts.append([geneToAdd, 1])
#
#
# print(geneCounts)
# with open(outFileName) as outFile:
#     for geneCount in geneCounts:
#         outFile.write(geneCount[0] + "\t" + str(geneCount[1]) + "\n")
