# gb contig files combined into a single fasta file. assumes that contigs are named like bioproject files "ADSFADF0100000xxx.gb"
# makes an intermediate file called temp.fasta so it will write over that file if it is in outputPath
from functions import *
pathToContigFiles = "./DownloadingFilesFromNCBI/AllMastitisFromSSB/"
outputPath = "./DownloadingFilesFromNCBI/AllMastitisFromSSB/CombinedFastas/"
if not os.path.exists(outputPath):
    print("creating output path")
    os.mkdir(outputPath)

with open(outputPath + "temp.fasta", "w") as tempOutFile:
    for filePath in glob(pathToContigFiles + "*.fasta"):
        numContigs = 0
        contig = ""
        with open(filePath) as file:
            for line in file:
                if line[0] == ">":
                    numContigs += 1
                else:
                    contig += line.strip()
        if numContigs > 1: # since it is a file of a single contig
            raise Exception("file has more than one contig" + str(numContigs))

        # look for locus name
        with open(filePath) as file:
            for line in file:
                cols = line.split()
                if len(cols) < 1:
                    continue
                if line[0] == ">":
                    tempOutFile.write(cols[0] + " whole genome shotgun sequence" + "\n") # write out the locus name
                    break
        tempOutFile.write(contig + "\n")

separateGenomesFromSingleFastaFile(outputPath + "temp.fasta")
os.remove(outputPath + "temp.fasta")