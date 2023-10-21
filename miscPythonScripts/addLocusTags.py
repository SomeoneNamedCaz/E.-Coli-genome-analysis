inFileName = "/Users/cazcullimore/Downloads/M12 closed genome.gb"
outFileName = "/Users/cazcullimore/Downloads/M12LocusTagged.gb"
indexOfDNARegionName = 0

with open(inFileName, "r") as fileToRead, open(outFileName, "w") as outFile:
    geneIndex = 0
    for line in fileToRead:
        cols = line.split()
        if cols[indexOfDNARegionName] == "gene" or cols[indexOfDNARegionName] == "CDS":
            line += "                     /locus_tag=b" + str(geneIndex) + "\n"
            geneIndex += 1
        outFile.write(line)
