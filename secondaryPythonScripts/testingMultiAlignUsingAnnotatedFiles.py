from functions import *
"""take all qor genes and find the distribution of nucs in each metadata catagory"""

geneName = "qorB"
pathToGBFiles = "/Users/cazcullimore/Documents/ericksonLabCode/filesToTestMultiAlign/gbks/*/*.gbk"
metadataFilePath = "/Users/cazcullimore/Documents/ericksonLabCode/metaDataForMetaCatsWithExtraMastitis.tsv"
fileNameToMetadataCategory = {}
with open(metadataFilePath) as metadataFile: # process metadata
    isFirstLine = True
    for line in metadataFile:
        line = line.strip()
        cols = line.split("\t")
        if isFirstLine:
            isFirstLine = False
        else:
            fileNameToMetadataCategory[cols[0].split(".")[0]] = cols[1:]
# outFile = "qorB."
snpIndex = 150
qorBs = {"pathogen": {"A":0,"T":0,"C":0,"G":0}, "commensal": {"A":0,"T":0,"C":0,"G":0}}

for path in glob(pathToGBFiles):
    try:
        fileName = path.split("/")[-1]
        contigs = getContigs(path)
        try:
            nuc = getGenesOnContigs(path, contigs)[geneName][1][snpIndex]
            print(nuc)
        except KeyError:
            print(fileName, "didn't have the gene", geneName)
        metadataCatagoryForFile = fileNameToMetadataCategory["scaffold_"+ fileName.split(".")[0]][1]
        try:
            qorBs[metadataCatagoryForFile][nuc] += 1
        except:
            print(nuc, "not in metadata")
    except:
        print(path, "unknown error")

print(qorBs)