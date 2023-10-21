import sys

snpGenomePath = sys.argv[1]#"./k-12RefGenomeAnalysisFiles/snpGenomeFiles/allSnps.afa"
mastitisMetaDataFromEnterobase = sys.argv[2] #"./allStrainsFromEntero.tsv"

# outputs
outputMetadataFileName = "mastitisSeverityMetadata.tsv"
outputSnpGenomeFileName = "allMastitisSeveritySnps.afa"

genomeNameToSequence = {}

with open(snpGenomePath) as file:
    currGenomeName = ""
    for line in file:
        line = line.strip()
        if line == "":
            continue

        if line[0] == ">":
            currGenomeName = line[1:]
        else:
            genomeNameToSequence[currGenomeName] = line


with open(mastitisMetaDataFromEnterobase) as inFile, open(outputMetadataFileName, "w") as metaOutFile, open(outputSnpGenomeFileName, "w") as snpOutFile:
    metaOutFile.write("strainName\tseverity\n")
    for line in inFile:
        line = line.strip()
        cols = line.split("\t")
        if line == "" or cols[0] == "Uberstrain" or cols[2] == "":
            continue
        barCode = cols[-1]
        severity = int(cols[2])
        if severity == 3:
            severity = "mild"
        elif severity == 1:
            severity = "severe"

        genomeName = "scaffold_"+ barCode + ".result.fasta.vcf"
        if genomeName in genomeNameToSequence.keys():
            metaOutFile.write("\t".join([genomeName, severity])  + "\n")
            snpOutFile.write(">" + genomeName + "\n" + genomeNameToSequence[genomeName] + "\n")
    metaOutFile.write("\n")