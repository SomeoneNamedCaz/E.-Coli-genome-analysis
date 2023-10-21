import os
from glob import glob
import re

fileFolderPath = "./assemblies"
tempFileName = "temp.txt"
for fileName in glob(os.path.join(fileFolderPath, '*.fasta')):
    with open(fileName) as assembly, open(tempFileName, "w") as tempFile:
        for line in assembly:
            if line.strip() != "":
                tempFile.write(line.strip() + "\n")
            else:
                print("removed line")
    with open(fileName, "w") as assembly, open(tempFileName) as tempFile:
        for line in tempFile:
            assembly.write(line)