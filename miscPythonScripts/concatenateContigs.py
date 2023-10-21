from functions import *

locationOfInFiles = "./allMastitisAssemblies/"
outPath = "./concatenatedMastitisAssemblies/"
for filePath in glob(locationOfInFiles + "*.fasta"):
    fileName = filePath.split("/")[-1]
    print(fileName)
    with open(filePath) as inFile, open(outPath + fileName, "w") as outFile:
        firstLine = True
        for line in inFile:
            if line[0] != ">" or firstLine:
                outFile.write(line)
                firstLine = False