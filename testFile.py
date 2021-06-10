from functions import *
pathBeginning = "./AllAssemblies/AllCommensalBovineAssemblies"

if (not os.path.exists(pathBeginning + "/ragtagOutputs/longestScaffoldFiles/")):
    os.mkdir(pathBeginning + "/ragtagOutputs/longestScaffoldFiles/")

for fileName in glob(pathBeginning + "/ragtagOutputs/*/ragtag.scaffolds.fasta"):
    uniqueName = fileName.split("/")[-2]
    with open(fileName) as file, open(pathBeginning + "/ragtagOutputs/longestScaffoldFiles/scaffold_" + uniqueName, "w") as outFile:
        for line in file:
            outFile.write(line)
            if line[0] != ">":
                # longestContig = line
                print(len(line))
                break

# for fileName in glob("./AllAssemblies/GCF_003018575.1_ASM301857v1_genomicCopy.fna"):
#     uniqueName = fileName.split("/")[-1]
#     with open(fileName) as file, open("./AllAssemblies/refSeqGCF_003018575.1_ASM301857JustChromosome.fna", "w") as outFile:
#         firstLine = True
#         for line in file:
#             outFile.write(line)
#             if line[0] == ">" and not firstLine:
#                 # longestContig = line
#                 # print(len(line))
#                 break
#             firstLine = False
