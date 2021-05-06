# gb contig files combined into a single fasta file. assumes that contigs are named like bioproject files "ADSFADF0100000xxx.gb"
# makes an intermediate file called temp.fasta so it will write over that file if it is in outputPath
from functions import *
pathToContigFiles = "/Users/cazcullimore/Documents/ericksonLabCode/DownloadingFilesFromNCBI/DownloadedFromSSB/"
outputPath = "/Users/cazcullimore/Documents/ericksonLabCode/DownloadingFilesFromNCBI/DownloadedFromSSB/CombinedFastas/"
if not os.path.exists(outputPath):
    print("creating output path")
    os.mkdir(outputPath)

with open(outputPath + "temp.fasta", "w") as tempOutFile:
    for filePath in glob(pathToContigFiles + "*.gb"):
        contig = getContigs(filePath)
        if len(contig) > 1: # since it is a list of contigs at this point
            raise Exception("file has more than one contig")
        contig = contig[0]
        # look for locus name
        with open(filePath) as file:
            for line in file:
                cols = line.split()
                if len(cols) < 1:
                    continue
                if cols[0] == "LOCUS":
                    tempOutFile.write(">" + cols[1] + " whole genome shotgun sequence" + "\n") # write out the locus name
                    break
        tempOutFile.write(contig + "\n")

separateGenomesFromSingleFastaFile(outputPath + "temp.fasta")
os.remove(outputPath + "temp.fasta")