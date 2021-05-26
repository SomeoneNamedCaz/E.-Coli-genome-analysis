from functions import *
# with open("./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/AllFilesInOneTry2.fasta" , "w") as outFile:
#     for fileName in glob("./allMastitisAssemblies/ragtagOutputs/*/ragtag.scaffolds.fasta"):
#         uniqueName = fileName.split("/")[-2]
#         with open(fileName) as file:
#             for line in file:
#                 if line[0] != ">":
#                     # longestContig = line
#                     # print(len(line))
#                     outFile.write(line)
#                     break
#                 else:
#                     outFile.write(">" + uniqueName + "\n")

for fileName in glob("./allMastitisAssemblies/GCF_003018575.1_ASM301857v1_genomicCopy.fna"):
    uniqueName = fileName.split("/")[-1]
    with open(fileName) as file, open("./refSeqGCF_003018575.1_ASM301857JustChromosome.fna", "w") as outFile:
        firstLine = True
        for line in file:
            outFile.write(line)
            if line[0] == ">" and not firstLine:
                # longestContig = line
                # print(len(line))
                break
            firstLine = False
