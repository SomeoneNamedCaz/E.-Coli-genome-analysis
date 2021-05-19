from functions import *

for fileName in glob("./allMastitisAssemblies/ragtagOutputs/*/ragtag.scaffolds.fasta"):
    uniqueName = fileName.split("/")[-2]
    with open(fileName) as file, open("./allMastitisAssemblies/ragtagOutputs/longestScaffoldFiles/scaffold_" + uniqueName, "w") as outFile:
        for line in file:
            outFile.write(line)
            if line[0] != ">":
                # longestContig = line
                # print(len(line))
                break
